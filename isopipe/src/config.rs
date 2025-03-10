use log::{error, info};
use serde::Deserialize;

use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::{Command, ExitStatus};
use std::thread;

use crate::cli::StepArgs;

const ISOTOOLS: &str = "isotools";

/// A struct representing a configuration file.
///
/// # Fields
///
/// * `metadata` - A HashMap containing metadata key-value pairs.
/// * `packages` - A HashMap containing package names.
/// * `steps` - A Vec containing PipelineStep enums.
///
/// # Example
///
/// ``` toml
/// [metadata]
/// key = "value"
///
/// [packages]
/// ccs = "5.0.0"
///
/// [[steps]]
/// step = "ccs"
/// ```
///
/// ``` rust, no_run
/// let config = Config {
///   metadata: HashMap::new(),
///   packages: HashMap::new(),
///   steps: vec![PipelineStep::Ccs],
/// };
/// ```
#[derive(Deserialize, Debug, Clone)]
pub struct Config {
    metadata: HashMap<String, String>,
    packages: HashMap<String, String>,
    #[serde(default, deserialize_with = "deserialize_steps")]
    steps: Vec<PipelineStep>,
}

impl Config {
    /// Read a configuration file and return a Config struct.
    ///
    /// # Arguments
    ///
    /// * `config` - A PathBuf containing the path to the configuration file.
    ///
    /// # Returns
    ///
    /// A Result containing a Config struct or an error.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let config = Config::read(PathBuf::from("config.toml"));
    /// ```
    pub fn read(config: PathBuf) -> Result<Self, Box<dyn std::error::Error>> {
        let mut file = File::open(config)?;
        let mut contents = String::new();
        file.read_to_string(&mut contents)?;

        let config: Config = toml::from_str(&contents)?;

        Ok(config)
    }

    /// Create a new Config struct.
    ///
    /// # Returns
    ///
    /// A Config struct.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let config = Config::new();
    /// ```
    pub fn new() -> Self {
        Self {
            metadata: HashMap::new(),
            packages: HashMap::new(),
            steps: Vec::new(),
        }
    }

    /// Config metadata getter.
    ///
    /// # Returns
    ///
    /// A reference to the metadata HashMap.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let config = Config::new();
    /// let metadata = config.metadata();
    /// ```
    pub fn metadata(&self) -> &HashMap<String, String> {
        &self.metadata
    }

    /// Config packages getter.
    ///
    /// # Returns
    ///
    /// A reference to the packages Vec.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let config = Config::new();
    /// let packages = config.packages();
    /// ```
    pub fn packages(&self) -> &HashMap<String, String> {
        &self.packages
    }

    /// Config steps getter.
    ///
    /// # Returns
    ///
    /// A reference to the steps Vec.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let config = Config::new();
    /// let steps = config.steps();
    /// ```
    pub fn steps(&self) -> &Vec<PipelineStep> {
        &self.steps
    }

    /// Config metadata setter.
    ///
    /// # Arguments
    ///
    /// * `metadata` - A HashMap containing metadata key-value pairs.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let mut config = Config::new();
    /// let mut metadata = HashMap::new();
    ///
    /// metadata.insert("key".into(), "value".into());
    ///
    /// config.set_metadata(metadata);
    /// ```
    pub fn set_metadata(&mut self, metadata: HashMap<String, String>) {
        self.metadata = metadata;
    }

    /// Config packages setter.
    ///
    /// # Arguments
    ///
    /// * `packages` - A Vec containing package names.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let mut config = Config::new();
    /// let packages = vec!["ccs".into(), "lima".into()];
    ///
    /// config.set_packages(packages);
    /// ```
    pub fn set_packages(&mut self, packages: HashMap<String, String>) {
        self.packages = packages;
    }

    /// Config steps setter.
    ///
    /// # Arguments
    ///
    /// * `steps` - A Vec containing PipelineStep enums.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let mut config = Config::new();
    /// let steps = vec![PipelineStep::Ccs, PipelineStep::Lima];
    ///
    /// config.set_steps(steps);
    /// ```
    pub fn set_steps(&mut self, steps: Vec<PipelineStep>) {
        self.steps = steps;
    }

    /// Add a package to the Config.
    ///
    /// # Arguments
    ///
    /// * `package` - A String containing the package name.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let mut config = Config::new();
    ///
    /// config.add_package("ccs".into());
    /// ```
    pub fn add_package(&mut self, package: String, version: String) {
        self.packages.insert(package, version);
    }

    /// Update packages in the Config based on the updated steps
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let mut config = Config::new();
    ///
    /// config.update_package();
    /// ```
    pub fn update_packages(&mut self) {
        let steps = &self.steps;

        for pkg in self.packages.clone().keys() {
            if !steps.iter().any(|&s| {
                if &pkg == &"isoseq3" {
                    return PipelineStep::Refine == s || PipelineStep::Cluster == s;
                }

                if &pkg == &"isotools" {
                    return PipelineStep::FilterQuality == s;
                }

                if &pkg == &"pbcss" {
                    return PipelineStep::Ccs == s;
                }

                s == PipelineStep::from_str(&pkg)
                    .expect("ERROR: Could not parse step from package name!")
            }) {
                self.packages.remove(pkg);
            }
        }
    }

    /// Add a step to the Config.
    ///
    /// # Arguments
    ///
    /// * `step` - A PipelineStep enum.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let mut config = Config::new();
    ///
    /// config.add_step(PipelineStep::Ccs);
    /// ```
    pub fn add_step(&mut self, step: PipelineStep) {
        self.steps.push(step);
    }

    /// Remove a package from the Config.
    ///
    /// # Arguments
    ///
    /// * `package` - A String containing the package name.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let mut config = Config::new();
    ///
    /// config.remove_package("ccs".into());
    /// ```
    pub fn remove_package(&mut self, package: String) {
        self.packages.remove(&package);
    }

    /// Remove a step from the Config.
    ///
    /// # Arguments
    ///
    /// * `step` - A PipelineStep enum.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let mut config = Config::new();
    ///
    /// config.remove_step(PipelineStep::Ccs);
    /// ```
    pub fn remove_step(&mut self, step: PipelineStep) {
        self.steps.retain(|s| s != &step);
    }

    /// Load packages from the Config.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let config = Config::new();
    /// config.load().unwrap();
    /// ```
    pub fn load(&self) -> Result<(), Box<dyn std::error::Error>> {
        for (package, version) in &self.packages {
            if package == "isotools" {
                build_isotools().expect("ERROR: Could not build isotools!");
            } else {
                load_package(package.clone(), Some(version.clone()))?;
            }
        }

        Ok(())
    }

    /// In-place modification of the steps in the Config
    ///
    /// # Arguments
    ///
    /// * `args` - A StepArgs struct containing the arguments.
    ///
    /// # Returns
    ///
    /// A mutable reference to the Config struct.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let mut config = Config::new();
    /// let args = StepArgs {
    ///    config: PathBuf::from("config.toml"),
    ///  from: "0".to_string(),
    /// to: "6".to_string(),
    /// only: None,
    /// skip: "3",
    /// dry_run: false,
    /// verbose: false,
    /// quiet: false,
    /// };
    ///
    /// config.aware(args);
    ///
    /// assert_eq!(config.steps().len(), 6);
    /// ```
    pub fn aware(&mut self, args: StepArgs) -> &mut Self {
        let steps = args
            .abs_steps()
            .expect("ERROR: An error ocurred while materializing steps!");

        self.set_steps(steps);
        self.update_packages();

        self
    }
}

/// Load a package using the module system.
///
/// # Arguments
///
/// * `package` - A String containing the package name.
/// * `version` - An Option containing the package version.
///
/// # Returns
///
/// A Result containing a unit or an error.
///
/// # Example
///
/// ``` rust, no_run
/// load_package("ccs".into(), Some("5.0.0".into()));
/// ```
fn load_package(
    package: String,
    version: Option<String>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = format!("module load {}", package);
    if let Some(v) = version {
        cmd.push_str(&format!("/{}", v));
    }
    log::info!("Loading: {}...", cmd);

    let output = std::process::Command::new("bash")
        .arg("-c")
        .arg(&cmd)
        .output()
        .expect("ERROR: Failed to execute process");

    if output.status.success() {
        Ok(())
    } else {
        log::error!(
            "ERROR: failed to load package {}\n{}.\nINFO: Trying locally...",
            cmd,
            String::from_utf8_lossy(&output.stderr)
        );
        // INFO: trying now locally to see if it works
        let output = std::process::Command::new(package.clone())
            .arg("--help")
            .output()
            .expect("ERROR: Failed to execute process");

        if output.status.success() {
            log::info!(
                "Package {} found: {}",
                package,
                String::from_utf8_lossy(&output.stdout)
            );
            Ok(())
        } else {
            log::error!(
                "ERROR: failed to execute {}\n{}\nPlease install the package.",
                package,
                String::from_utf8_lossy(&output.stderr)
            );
            std::process::exit(1);
        }
    }
}

impl Default for Config {
    fn default() -> Self {
        Self::new()
    }
}

/// An enum representing pipeline steps.
///
/// # Example
///
/// ``` rust, no_run
/// let step = PipelineStep::Ccs;
/// ```
#[derive(Deserialize, Debug, Clone, Copy, PartialEq, PartialOrd, Eq, Ord)]
pub enum PipelineStep {
    Ccs,
    Lima,
    Refine,
    Cluster,
    Minimap,
    FilterQuality,
    LoadGenome,
}

impl PipelineStep {
    /// Create a PipelineStep enum from a string.
    ///
    /// # Arguments
    ///
    /// * `s` - A string containing the pipeline step.
    ///
    /// # Returns
    ///
    /// A Result containing a PipelineStep enum or an error.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let step = PipelineStep::from_str("ccs");
    ///
    /// assert_eq!(step, Ok(PipelineStep::Ccs));
    /// ```
    pub fn from_str(s: &str) -> Result<Self, String> {
        match s {
            "ccs" => Ok(Self::Ccs),
            "lima" => Ok(Self::Lima),
            "refine" => Ok(Self::Refine),
            "cluster" => Ok(Self::Cluster),
            "minimap3" => Ok(Self::Minimap),
            "filter-quality" => Ok(Self::FilterQuality),
            "load-genome" => Ok(Self::LoadGenome),
            _ => Err(format!("ERROR: Invalid pipeline step: {}", s)),
        }
    }

    /// Create a PipelineStep enum from an integer.
    ///
    /// # Arguments
    ///
    /// * `i` - An integer containing the pipeline step.
    ///
    /// # Returns
    ///
    /// A Result containing a PipelineStep enum or an error.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let step = PipelineStep::from_int(0);
    ///
    /// assert_eq!(step, Ok(PipelineStep::Ccs));
    /// ```
    pub fn from_int(i: usize) -> Result<Self, String> {
        match i {
            0 => Ok(Self::Ccs),
            1 => Ok(Self::Lima),
            2 => Ok(Self::Refine),
            3 => Ok(Self::Cluster),
            4 => Ok(Self::Minimap),
            5 => Ok(Self::FilterQuality),
            6 => Ok(Self::LoadGenome),
            _ => Err(format!("ERROR: Invalid pipeline step: {}", i)),
        }
    }

    /// Convert a PipelineStep enum to a string.
    ///
    /// # Returns
    ///
    /// A string containing the pipeline step.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let step = PipelineStep::Ccs;
    /// let s = step.to_str();
    ///
    /// assert_eq!(s, "ccs");
    /// ```
    pub fn to_str(&self) -> String {
        match self {
            Self::Ccs => "ccs".into(),
            Self::Lima => "lima".into(),
            Self::Refine => "refine".into(),
            Self::Cluster => "cluster".into(),
            Self::Minimap => "minimap3".into(),
            Self::FilterQuality => "filter-quality".into(),
            Self::LoadGenome => "load-genome".into(),
        }
    }

    /// Convert a PipelineStep enum to an integer.
    ///
    /// # Returns
    ///
    /// An integer containing the pipeline step.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let step = PipelineStep::Ccs;
    /// let i = step.to_int();
    ///
    /// assert_eq!(i, 0);
    /// ```
    pub fn to_int(&self) -> usize {
        match self {
            Self::Ccs => 0,
            Self::Lima => 1,
            Self::Refine => 2,
            Self::Cluster => 3,
            Self::Minimap => 4,
            Self::FilterQuality => 5,
            Self::LoadGenome => 6,
        }
    }

    /// Convert a Vec of strings to a Vec of PipelineStep enums.
    ///
    /// # Arguments
    ///
    /// * `v` - A Vec of strings.
    ///
    /// # Returns
    ///
    /// A Result containing a Vec of PipelineStep enums or an error.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let v = vec!["ccs".into(), "lima".into()];
    /// let steps = PipelineStep::from_vec_str(v);
    ///
    /// assert_eq!(steps, Ok(vec![PipelineStep::Ccs, PipelineStep::Lima]));
    /// ```
    pub fn from_vec_str(v: Vec<String>) -> Result<Vec<Self>, String> {
        v.iter().map(|s| Self::from_str(s)).collect()
    }

    /// Convert a Vec of integers to a Vec of PipelineStep enums.
    ///
    /// # Arguments
    ///
    /// * `v` - A Vec of integers.
    ///
    /// # Returns
    ///
    /// A Result containing a Vec of PipelineStep enums or an error.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let v = vec![0, 1];
    /// let steps = PipelineStep::from_vec_int(v);
    ///
    /// assert_eq!(steps, Ok(vec![PipelineStep::Ccs, PipelineStep::Lima]));
    /// ```
    pub fn from_vec_int(v: Vec<usize>) -> Result<Vec<Self>, String> {
        v.iter().map(|i| Self::from_int(*i)).collect()
    }
}

/// Deserialize a Vec of PipelineStep enums.
///
/// # Arguments
///
/// * `deserializer` - A serde Deserializer.
///
/// # Returns
///
/// A Result containing a Vec of PipelineStep enums or an error.
///
/// # Example
///
/// ``` rust, no_run
/// let steps = deserialize_steps("ccs,lima");
/// ```
///
/// ``` toml
/// steps = "ccs,lima"
/// ```
fn deserialize_steps<'de, D>(deserializer: D) -> Result<Vec<PipelineStep>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    let s: Option<String> = Option::deserialize(deserializer)?;
    Ok(s.map_or(vec![], |_| vec![]))
}

/// Run a command and return the exit status.
///
/// # Arguments
///
/// * `cmd` - A mutable reference to a Command.
///
/// # Returns
///
/// An ExitStatus.
///
/// # Example
///
/// ``` rust, no_run
/// let mut cmd = Command::new("ls");
/// let status = run_command(&mut cmd);
///
/// assert!(status.success());
/// ```
fn run_command(cmd: &mut Command) -> ExitStatus {
    match cmd.status() {
        Ok(status) => status,
        Err(e) => {
            error!("ERROR: Failed to execute process: {}", e);
            std::process::exit(1);
        }
    }
}

/// Builds the isotools submodule.
///
/// # Example
///
/// ``` rust, no_run
/// build_isotools();
///
/// assert!(run_command(Command::new("isotools").arg("--help")).success());
/// ```
fn build_isotools() -> Result<(), Box<dyn std::error::Error>> {
    let isotools = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("ERROR: Could not find parent directory")
        .join(ISOTOOLS);

    // INFO: check if isotools is already installed
    if run_command(Command::new("isotools").arg("--help")).success() {
        info!("isotools already installed!");
        return Ok(());
    }

    // INFO: onitialize and update the submodule in one step
    if !run_command(Command::new("git").args(["submodule", "update", "--init"])).success() {
        error!("ERROR: Failed to initialize/update isotools submodule.");
        std::process::exit(1);
    }

    info!("Submodule initialized and updated successfully.");

    // INFO: build and install isotools concurrently
    let build_path = isotools.clone();
    let install_path = isotools.join("entry");

    let build_thread = thread::spawn(move || cargo_build(&build_path));
    let install_thread = thread::spawn(move || cargo_install(&install_path));

    build_thread.join().expect("ERROR: Build thread panicked.");
    install_thread
        .join()
        .expect("ERROR: Install thread panicked.");

    // INFO: final check if isotools is installed
    if run_command(Command::new("isotools").arg("--help")).success() {
        info!("isotools successfully installed.");
    } else {
        error!("ERROR: Failed to install isotools.");
        std::process::exit(1);
    }

    Ok(())
}

/// Builds a cargo project.
///
/// # Arguments
///
/// * `path` - A PathBuf containing the path to the cargo project.
///
/// # Example
///
/// ``` rust, no_run
/// cargo_build(PathBuf::from("isotools"));
///
/// assert!(run_command(Command::new("isotools").arg("--help")).success());
/// ```
fn cargo_build(path: &Path) {
    if !run_command(Command::new("cargo").arg("build").current_dir(path)).success() {
        error!("ERROR: Failed to build cargo project!");
        std::process::exit(1);
    }

    info!(
        "{}",
        format!("Cargo project from {} build complete.", path.display())
    );
}

/// Installs a cargo project.
///
/// # Arguments
///
/// * `path` - A PathBuf containing the path to the cargo project.
///
/// # Example
///
/// ``` rust, no_run
/// cargo_install(PathBuf::from("isotools"));
///
/// assert!(run_command(Command::new("isotools").arg("--help")).success());
/// ```
fn cargo_install(path: &Path) {
    if !run_command(Command::new("cargo").args(["install", "--path"]).arg(path)).success() {
        error!("ERROR: Failed to install cargo project!");
        std::process::exit(1);
    }

    info!(
        "{}",
        format!(
            "Cargo project from {} installation complete.",
            path.display()
        )
    );
}
