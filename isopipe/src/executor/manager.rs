use std::io::Write;
use std::str::FromStr;

use crate::{
    config::{Config, PipelineStep},
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
    /// executor.execute();
    /// ```
    pub fn execute(&self, config: &Config, step: &PipelineStep) {
        let jobs = write_jobs(self.jobs.clone());

        match self.manager {
            ParallelManager::Nextflow => {
                // INFO: 'nextflow run <pipeline> -c <config> -j <jobs>'
                let cmd = format!(
                    "nextflow run --jobs {} --mem {}",
                    jobs,
                    config.get_param(*step, "memory").to_int()
                );

                std::process::Command::new("sh")
                    .arg("-c")
                    .arg(cmd)
                    .output()
                    .expect("ERROR: Failed to execute command");
            }
            ParallelManager::Para => {
                // INFO: 'para make <step> <jobs> -q <queue> -memoryMb <memory>'
                let cmd = format!(
                    "para make {} {} -q {} -memoryMb {}",
                    step,
                    jobs,
                    config
                        .global
                        .get("short_queue")
                        .expect("ERROR: No short queue found"),
                    config.get_param(*step, "memory").to_int(),
                );

                std::process::Command::new("sh")
                    .arg("-c")
                    .arg(cmd)
                    .output()
                    .expect("ERROR: Failed to execute command");
            }
            ParallelManager::Snakemake => {
                todo!()
            }
            ParallelManager::Local => {
                todo!()
            }
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
    pub fn init(manager: &str) -> Self {
        check_manager(manager);

        match manager.to_lowercase().as_str() {
            "nextflow" => ParallelManager::Nextflow,
            "para" => ParallelManager::Para,
            "snakemake" => ParallelManager::Snakemake,
            "local" => ParallelManager::Local,
            _ => panic!("ERROR: Unknown parallel manager: {}", manager),
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
fn write_jobs(jobs: Vec<Job>) -> String {
    let filename = String::from("jobs");

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
fn check_manager(manager: &str) {
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
