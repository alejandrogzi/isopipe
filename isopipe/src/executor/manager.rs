use std::io::Write;
use std::path::PathBuf;
use std::str::FromStr;

use crate::{
    config::{Config, PipelineStep},
    consts::*,
    executor::job::Job,
};

#[derive(Debug, Clone)]
pub struct ParallelExecutor {
    /// Command to run the parallel manager
    pub manager: ParallelManager,
    /// List of jobs to run
    pub jobs: Vec<Job>,
    /// List of arguments to pass to the parallel manager
    pub args: Vec<String>,
}

impl ParallelExecutor {
    /// Create a new instance of ParallelExecutor
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use isopipe::executor::ParallelExecutor;
    ///
    /// let executor = ParallelExecutor::new(ParallelManager::Nextflow);
    /// assert_eq!(executor.manager, "nextflow");
    /// ```
    pub fn new(manager: ParallelManager) -> Self {
        Self {
            manager,
            jobs: Vec::new(),
            args: Vec::new(),
        }
    }

    /// Add a job to the executor
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use isopipe::executor::ParallelExecutor;
    ///
    /// let executor = ParallelExecutor::new(ParallelManager::Nextflow);
    /// executor.add_job("job1");
    /// ```
    pub fn add_job(&mut self, job: Job) {
        self.jobs.push(job);
    }

    /// Add a list of jobs to the executor
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use isopipe::executor::ParallelExecutor;
    ///
    /// let executor = ParallelExecutor::new(ParallelManager::Nextflow);
    /// executor.add_jobs(vec!["job1", "job2"]);
    /// ```
    pub fn add_jobs(&mut self, jobs: Vec<Job>) -> &mut Self {
        self.jobs.extend(jobs);

        self
    }

    /// Add a list of arguments to the executor
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use isopipe::executor::ParallelExecutor;
    ///
    /// let executor = ParallelExecutor::new(ParallelManager::Nextflow);
    /// executor.add_args(vec!["--arg1".to_string(), "--arg2".to_string()]);
    /// ```
    pub fn add_args(&mut self, args: Vec<String>) {
        self.args.extend(args);
    }

    /// Add an argument to the executor
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use isopipe::executor::ParallelExecutor;
    ///
    /// let executor = ParallelExecutor::new(ParallelManager::Nextflow);
    /// executor.add_arg("--arg1".to_string());
    /// ```
    pub fn add_arg(&mut self, arg: String) {
        self.args.push(arg);
    }

    /// Execute the parallel manager with the jobs and arguments
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use isopipe::executor::ParallelExecutor;
    ///
    /// let executor = ParallelExecutor::new(ParallelManager::Nextflow);
    /// executor.add_job("job1");
    /// executor.add_arg("--arg1".to_string());
    ///
    /// assert_eq!(executor.manager, "nextflow");
    /// assert_eq!(executor.jobs.len(), 1);
    /// assert_eq!(executor.args.len(), 1);
    ///
    /// executor.execute();
    /// ```
    pub fn execute(&mut self, config: &Config, step: &PipelineStep, global_output_dir: PathBuf) {
        let jobs = write_jobs(self.jobs.clone(), global_output_dir.clone());
        let package = config.get_package_from_step(step);

        let memory = config
            .get_param(*step, "memory")
            .unwrap_or(
                config
                    .get_global_param(DEFAULT_MEMORY)
                    .expect("ERROR: No default memory found in global parameters!"),
            )
            .to_int();

        let threads = config
            .get_param(*step, "num-threads")
            .unwrap_or(
                // INFO: covering minimap -t flag -> else, roll up to default
                config.get_param(*step, "t").unwrap_or(
                    config
                        .get_global_param(DEFAULT_THREADS)
                        .expect("ERROR: No default threads found in global parameters!"),
                ),
            )
            .to_int();

        match self.manager {
            ParallelManager::Nextflow => {
                // INFO: 'nextflow run <pipeline> -j <jobs>'
                let runner = __get_assets_dir().join(NF_RUNNER);

                let cmd = format!(
                    "module load {} && nextflow run {} --jobs {} --mem {} --threads {}",
                    package,
                    runner.display(),
                    jobs.display(),
                    memory,
                    threads,
                );

                std::process::Command::new("sh")
                    .arg("-c")
                    .arg(cmd)
                    .output()
                    .expect("ERROR: Failed to execute command");
            }
            ParallelManager::Para => {
                // INFO: 'para make <step> <jobs> -q <queue> -memoryMb <memory>'
                self.__para(
                    config,
                    &step.to_unique_str(),
                    &jobs,
                    threads as u32,
                    memory as u32,
                    package,
                );
            }
            ParallelManager::Snakemake => {
                todo!()
            }
            ParallelManager::Local => {
                todo!()
            }
        }

        self.reset(global_output_dir);
    }

    /// Reset the executor by clearing the jobs and arguments
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use isopipe::executor::ParallelExecutor;
    ///
    /// let executor = ParallelExecutor::new(ParallelManager::Nextflow);
    /// executor.add_job("job1");
    /// executor.add_arg("--arg1".to_string());
    ///
    /// assert_eq!(executor.manager, "nextflow");
    /// assert_eq!(executor.jobs.len(), 1);
    ///
    /// executor.reset();
    ///
    /// assert_eq!(executor.jobs.len(), 0);
    /// assert_eq!(executor.args.len(), 0);
    /// ```
    pub fn reset(&mut self, global_output_dir: PathBuf) {
        self.jobs.clear();
        self.args.clear();

        // INFO: remove jobs file
        let filename = global_output_dir.join("jobs");
        std::fs::remove_file(&filename).expect("ERROR: Failed to remove job file");
    }

    /// Channels errors while using para as executor
    ///
    /// # Arguments
    /// * `step` - The pipeline step that caused the error
    ///
    /// # Examples
    /// ```
    /// let step = PipelineStep::new("test");
    /// let executor = ParallelExecutor::new(ParallelManager::Para);
    /// executor.__channel_error(&step);
    /// ```
    pub fn __channel_error(&self, step: &PipelineStep, run_id: String) {
        match self.manager {
            ParallelManager::Para => {
                if std::env::current_dir()
                    .expect("ERROR: Failed to get current directory!")
                    .join(".para")
                    .exists()
                {
                    let dir = std::env::current_dir()
                        .expect("ERROR: Failed to get current directory")
                        .join(".para")
                        .join(format!("{}_{}", step, run_id))
                        .join("1");

                    log::info!(
                        "INFO: Checking for crashed processes in directory: {}...",
                        dir.display()
                    );

                    // INFO: looping through dir until find .crashed
                    for entry in std::fs::read_dir(&dir).expect("ERROR: Failed to read directory") {
                        let entry = entry.expect(&format!(
                            "ERROR: Failed to read directory entry -> {}",
                            dir.display(),
                        ));
                        if entry
                            .file_name()
                            .to_str()
                            .expect(&format!("ERROR: could not get filename -> {:?}", entry))
                            .ends_with(".crashed")
                        {
                            let error = std::fs::read_to_string(entry.path()).expect(&format!(
                                "ERROR: Failed to read error file -> {}",
                                entry.path().display()
                            ));
                            log::error!("ERROR: {}", error);

                            return;
                        }
                    }
                } else {
                    log::info!("INFO: starting clean run, no .para dir found!")
                }
            }
            _ => {}
        }
    }

    /// Send a custom set of jobs to the executor
    ///
    /// # Arguments
    ///
    /// * `config` - The configuration for the executor
    /// * `step` - The step to be executed
    /// * `dir` - The directory containing the jobs
    /// * `threads` - The number of threads to use
    /// * `memory` - The amount of memory to use
    /// * `package` - The package to be executed
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let config = Config::default();
    /// let step = "step1";
    /// let dir = PathBuf::from("/path/to/jobs");
    /// let threads = 4;
    /// let memory = 1024;
    /// let package = "package1";
    ///
    /// executor.and_send(&config, step, dir, threads, memory, package);
    /// ```
    pub fn and_send(
        &mut self,
        config: &Config,
        step: &str,
        dir: PathBuf,
        threads: u32,
        memory: u32,
        package: String,
    ) {
        match self.manager {
            ParallelManager::Para => {
                let jobs = write_jobs(self.jobs.clone(), dir);
                self.__para(config, step, &jobs, threads, memory, package);
            }
            _ => {
                todo!()
            }
        }
    }

    /// Send a custom set of jobs to para
    ///
    /// # Arguments
    /// * `config` - The configuration for the executor
    /// * `step` - The step name
    /// * `jobs` - The path to the jobs file
    /// * `threads` - The number of threads to use
    /// * `memory` - The amount of memory to use in MB
    /// * `package` - The package to use
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let config = Config::default();
    /// let step = "step1";
    /// let dir = PathBuf::from("/path/to/jobs");
    /// let threads = 4;
    /// let memory = 1024;
    /// let package = "package1";
    ///
    /// let mut manager = ExecutorManager::new(config);
    /// manager.__para(&config, step, &dir, threads, memory, package);
    /// ```
    pub fn __para(
        &mut self,
        config: &Config,
        step: &str,
        jobs: &PathBuf,
        threads: u32,
        memory: u32,
        package: String,
    ) {
        let run_id = config.get_run_id();
        let step_code = format!("{}_{}", step, run_id);

        let cmd = format!(
            "module load {} && para make {} {} -q {} -memoryMb {} -numCores {}",
            package,
            step_code,
            jobs.display(),
            config
                .global
                .get(SHORT_QUEUE)
                .expect("ERROR: No short queue found"),
            memory * 1024, // WARN: Memory is in MB
            threads,
        );

        log::info!("INFO: Executing command: {}", cmd);

        let output = std::process::Command::new("sh")
            .arg("-c")
            .arg(cmd)
            .output()
            .expect("ERROR: Failed to execute command");

        if !output.status.success() {
            log::error!(
                "ERROR: Failed to execute command: {}",
                String::from_utf8_lossy(&output.stderr)
            );

            if let Ok(step) = PipelineStep::from_str(step) {
                self.__channel_error(&step, run_id);
            }

            std::process::exit(1);
        } else {
            log::info!(
                "INFO: Command executed successfully: {}",
                String::from_utf8_lossy(&output.stdout)
            );
        }
    }
}

#[derive(Debug, Clone)]
pub enum ParallelManager {
    /// Parallel manager for Nextflow
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use isopipe::executor::ParallelManager;
    ///
    /// let manager = ParallelManager::Nextflow::new();
    /// assert_eq!(manager, "nextflow");
    /// ```
    Nextflow,

    /// Parallel manager for para
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use isopipe::executor::ParallelManager;
    ///
    /// let manager = ParallelManager::Para::new();
    /// assert_eq!(manager, "para");
    /// ```
    Para,

    /// Parallel manager for Snakemake
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use isopipe::executor::ParallelManager;
    ///
    /// let manager = ParallelManager::Snakemake::new();
    /// assert_eq!(manager, "snakemake");
    /// ```
    Snakemake,

    /// Local parallel manager
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use isopipe::executor::ParallelManager;
    ///
    /// let manager = ParallelManager::Local::new();
    /// assert_eq!(manager, "local");
    /// ```
    Local,
}

impl FromStr for ParallelManager {
    type Err = String;

    /// Convert a string to a ParallelManager
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use isopipe::executor::ParallelManager;
    ///
    /// let manager = ParallelManager::from_str("nextflow").unwrap();
    /// assert_eq!(manager, "nextflow");
    /// ```
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "nextflow" => Ok(ParallelManager::Nextflow),
            "para" => Ok(ParallelManager::Para),
            "snakemake" => Ok(ParallelManager::Snakemake),
            "local" => Ok(ParallelManager::Local),
            _ => Err(format!("ERROR: Unknown parallel manager: {}", s)),
        }
    }
}

impl From<&str> for ParallelManager {
    /// Convert a string to a ParallelManager
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use isopipe::executor::ParallelManager;
    ///
    /// let manager = ParallelManager::from_str("nextflow").unwrap();
    /// assert_eq!(manager, "nextflow");
    /// ```
    fn from(value: &str) -> Self {
        match value {
            "nextflow" => ParallelManager::Nextflow,
            "para" => ParallelManager::Para,
            "snakemake" => ParallelManager::Snakemake,
            "local" => ParallelManager::Local,
            _ => panic!("ERROR: Unknown parallel manager: {}", value),
        }
    }
}

impl std::fmt::Display for ParallelManager {
    /// Display the ParallelManager as a string
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use isopipe::executor::ParallelManager;
    ///
    /// let manager = ParallelManager::Nextflow;
    /// assert_eq!(format!("{}", manager), "nextflow");
    /// ```
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParallelManager::Nextflow => write!(f, "nextflow"),
            ParallelManager::Para => write!(f, "para"),
            ParallelManager::Snakemake => write!(f, "snakemake"),
            ParallelManager::Local => write!(f, "local"),
        }
    }
}

impl ParallelManager {
    /// Create a new instance of ParallelManager
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use isopipe::executor::ParallelManager;
    ///
    /// let manager = ParallelManager::new("nextflow");
    /// assert_eq!(manager.cmd, "nextflow");
    /// ```
    pub fn new(manager: &str) -> Self {
        __check_manager(manager);

        match manager.to_lowercase().as_str() {
            "nextflow" => ParallelManager::Nextflow,
            "para" => ParallelManager::Para,
            "snakemake" => ParallelManager::Snakemake,
            "local" => ParallelManager::Local,
            _ => panic!("ERROR: Unknown parallel manager: {}", manager),
        }
    }

    /// Initialize the parallel manager
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use isopipe::executor::ParallelManager;
    ///
    /// let manager = ParallelManager::new("nextflow");
    /// let executor = manager.init();
    ///
    /// assert_eq!(executor.manager, "nextflow");
    /// assert_eq!(executor.jobs.len(), 0);
    /// ```
    pub fn init(&self) -> ParallelExecutor {
        match self {
            ParallelManager::Nextflow => {
                // INFO: 'nextflow run <pipeline> -c <config> -j <jobs>'
                log::info!("INFO: Initializing nextflow...");
                self.as_executor()
            }
            ParallelManager::Para => {
                // INFO: 'para make <step> <jobs> -q <queue> -memoryMb <memory>'
                log::info!("INFO: Initializing para...");
                self.as_executor()
            }
            ParallelManager::Snakemake => {
                todo!()
            }
            ParallelManager::Local => {
                // INFO: 'sh -c <command>'
                log::info!("INFO: Initializing local...");
                self.as_executor()
            }
        }
    }

    /// Create a new instance of ParallelExecutor from the ParallelManager
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use isopipe::executor::ParallelManager;
    ///
    /// let manager = ParallelManager::new("nextflow");
    /// let executor = manager.as_executor();
    ///
    /// assert_eq!(executor.manager, "nextflow");
    /// assert_eq!(executor.jobs.len(), 0);
    /// ```
    pub fn as_executor(&self) -> ParallelExecutor {
        ParallelExecutor {
            manager: self.clone(),
            jobs: Vec::new(),
            args: Vec::new(),
        }
    }
}

/// Write the jobs to a file
///
/// # Example
///
/// ```rust, no_run
/// use isopipe::executor::write_jobs;
///
/// let jobs = vec![
///    Job::new("job1"),
///    Job::new("job2"),
/// ];
///
/// let filename = write_jobs(jobs);
///
/// assert_eq!(filename.to_str().unwrap(), "jobs");
/// ```
fn write_jobs(jobs: Vec<Job>, global_output_dir: PathBuf) -> PathBuf {
    let filename = global_output_dir.join("jobs");

    let mut file = std::fs::File::create(&filename).expect("ERROR: Failed to create job file");
    for job in jobs {
        let cmd = job.cmd;
        writeln!(file, "{}", cmd).expect("ERROR: Failed to write to job file");
    }

    filename
}

/// Check if the parallel manager is valid
///
/// # Example
///
/// ```rust, no_run
/// use isopipe::executor::check_manager;
///
/// let manager = "nextflow";
/// check_manager(manager);
///
/// assert_eq!(manager, "nextflow");
/// ```
fn __check_manager(manager: &str) {
    if !["nextflow", "para", "snakemake", "local"].contains(&manager) {
        panic!("ERROR: Unknown parallel manager: {}", manager);
    }

    if !manager.is_empty() {
        panic!("ERROR: Parallel manager cannot be empty");
    }

    std::process::Command::new(manager)
        .arg("--version")
        .output()
        .expect("ERROR: Failed to execute command");
}

/// Get the assets directory
///
/// # Example
///
/// ```rust, no_run
/// use isopipe::executor::get_assets_dir;
///
/// let assets_dir = get_assets_dir();
///
/// assert_eq!(assets_dir.to_str().unwrap(), "assets");
/// ```
pub fn __get_assets_dir() -> PathBuf {
    let assets = std::env::current_dir().expect("Failed to get executable path");
    return assets.join(ASSETS);
}

/// Clean the .para directory
///
/// # Arguments
/// * `dir` - The current directory where .para resides
///
/// # Example
///
/// ```rust, no_run
/// use isopipe::executor::__clean_para_dir;
///
/// let dir = std::env::current_dir().expect("Failed to get executable path");
/// let dir = __clean_para_dir(dir);
/// ```
#[allow(dead_code)]
pub fn __clean_para_dir(dir: PathBuf) {
    let para_dir = dir.join(".para");

    if para_dir.exists() {
        std::fs::remove_dir_all(&para_dir).expect("ERROR: Failed to remove para directory");
    }
}
