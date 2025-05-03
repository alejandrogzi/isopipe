use log::{error, info};
use serde::Deserialize;

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::{Command, ExitStatus};
use std::thread;

use crate::cli::StepArgs;
use crate::consts::*;

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
    pub metadata: HashMap<String, String>,
    pub packages: HashMap<String, String>,
    pub global: HashMap<String, ParamValue>,
    #[serde(default, deserialize_with = "deserialize_steps")]
    pub steps: Vec<PipelineStep>,
    #[serde(default, deserialize_with = "deserialize_to_hash")]
    pub params: HashMap<PipelineStep, StepParams>,
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
            global: HashMap::new(),
            steps: Vec::new(),
            params: HashMap::new(),
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
                if &pkg == &"isoseq" {
                    return PipelineStep::Refine == s || PipelineStep::Cluster == s;
                }

                if &pkg == &"isotools" {
                    return PipelineStep::Polya == s;
                }

                if &pkg == &"pbccs" || &pkg == &"pbindex" {
                    return PipelineStep::Ccs == s;
                }

                s == PipelineStep::from_str(&pkg)
                    .expect("ERROR: Could not parse step from package name!")
            }) {
                self.packages.remove(pkg);
            }

            if &pkg == &"pbcss" {
                // WARN: force pbindex if CCS is used
                if !self.packages.contains_key("pbindex") {
                    self.packages.insert("pbindex".into(), "1.7.0".into());
                }
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
    pub fn load(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        self.set_run_id();

        for (package, _) in &self.packages {
            if package == ISOTOOLS {
                build_isotools().expect("ERROR: Could not build isotools!");
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
    ///     config: PathBuf::from("config.toml"),
    ///     from: "0".to_string(),
    ///     to: "6".to_string(),
    ///     only: None,
    ///     skip: "3",
    ///     dry_run: false,
    ///     verbose: false,
    ///     quiet: false,
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
        self.update_params();

        self
    }

    /// Update the parameters in the Config based on the updated steps
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let mut config = Config::new();
    ///
    /// config.update_params();
    /// ```
    pub fn update_params(&mut self) {
        let steps = &self.steps;

        for step in self.params.clone().keys() {
            if !steps.iter().any(|&s| s == *step) {
                self.params.remove(step);
            }
        }
    }

    /// Config parameters getter.
    ///
    /// # Returns
    ///
    /// A reference to the parameters HashMap.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let config = Config::new();
    /// let params = config.params();
    /// ```
    pub fn params(&self) -> &HashMap<PipelineStep, StepParams> {
        &self.params
    }

    /// Config parameters setter.
    ///
    /// # Arguments
    ///
    /// * `params` - A HashMap containing step parameters.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let mut config = Config::new();
    /// let mut params = HashMap::new();
    ///
    /// params.insert(PipelineStep::Ccs, StepParams { values: HashMap::new() });
    ///
    /// config.set_params(params);
    ///
    /// assert_eq!(config.params().len(), 1);
    /// ```
    pub fn set_params(&mut self, params: HashMap<PipelineStep, StepParams>) {
        self.params = params;
    }

    /// Add a parameter to the Config.
    ///
    /// # Arguments
    ///
    /// * `step` - A PipelineStep enum.
    /// * `params` - A StepParams struct.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let mut config = Config::new();
    ///
    /// config.add_param(PipelineStep::Ccs, StepParams { values: HashMap::new() });
    /// ```
    pub fn add_param(&mut self, step: PipelineStep, params: StepParams) {
        self.params.insert(step, params);
    }

    /// Remove a parameter from the Config.
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
    /// config.remove_param(PipelineStep::Ccs, "min-rq".into());
    /// ```
    pub fn remove_param(&mut self, step: PipelineStep, key: &str) {
        let params = self
            .params
            .get_mut(&step)
            .expect("ERROR: Step not found in params!");

        params.values.remove(key);
    }

    /// Remove a step parameters from the Config.
    ///
    /// # Arguments
    ///
    /// * `step` - A PipelineStep enum.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let mut config = Config::new();
    /// config.remove_step_param(PipelineStep::Ccs);
    ///
    /// assert_eq!(config.params().get(&PipelineStep::Ccs), None);
    /// ```
    pub fn remove_step_param(&mut self, step: PipelineStep) {
        self.params.remove(&step);
    }

    /// Add a parameter value to a Config step.
    ///
    /// # Arguments
    ///
    /// * `step` - A PipelineStep enum.
    /// * `key` - A String containing the parameter key.
    /// * `value` - A ParamValue enum.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let mut config = Config::new();
    ///
    /// config.add_param_value(PipelineStep::Ccs, "min-rq".into(), ParamValue::Float(0.95));
    /// ```
    pub fn add_param_value(&mut self, step: PipelineStep, key: String, value: ParamValue) {
        let params = self
            .params
            .get_mut(&step)
            .expect("ERROR: Step not found in params!");

        params.values.insert(key, value);
    }

    /// Get a parameter value from a Config step.
    ///
    /// # Arguments
    ///
    /// * `step` - A PipelineStep enum.
    /// * `key` - A String containing the parameter key.
    ///
    /// # Returns
    ///
    /// A reference to the ParamValue enum.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let config = Config::new()
    ///   .add_param_value(PipelineStep::Ccs, "min-rq".into(), ParamValue::Float(0.95));
    /// let value = config.get_param(PipelineStep::Ccs, "min-rq");
    ///
    /// assert_eq!(value, ParamValue::Float(0.95));
    /// ```
    pub fn get_param(&self, step: PipelineStep, key: &str) -> Option<&ParamValue> {
        self.params
            .get(&step)
            .expect(format!("ERROR: Step {} not found in params!", step).as_str())
            .get(key)
    }

    /// Get a global parameter value from the Config.
    ///
    /// # Arguments
    ///
    /// * `key` - A String containing the parameter key.
    ///
    /// # Returns
    ///
    /// An Option containing the parameter value.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let config = Config::new();
    /// let value = config.get_global_param("min-rq");
    ///
    /// assert_eq!(value, None);
    /// ```
    pub fn get_global_param(&self, key: &str) -> Option<&ParamValue> {
        self.global.get(key)
    }

    /// Get global output directory from the Config
    /// with a timestamp appended.
    ///
    /// # Returns
    ///
    /// A PathBuf containing the global output directory.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let config = Config::new();
    /// let output = config.get_global_output();
    ///
    /// assert_eq!(output, PathBuf::from("output_20210901120000"));
    /// ```
    pub fn create_global_output_dir(&self) -> PathBuf {
        let rs = format!(
            "{}/{}_{}",
            self.global
                .get("global_output_dir")
                .expect("ERROR: output not found!")
                .to_path_buf()
                .display(),
            OUTPUT,
            chrono::Local::now().format("%Y%m%d%H%M")
        )
        .into();

        std::fs::create_dir_all(&rs).expect("ERROR: Could not create output directory!");

        rs
    }

    /// Get global data prefix from the Config.
    ///
    /// # Returns
    ///
    /// A String containing the global data prefix.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let config = Config::new();
    /// let data_prefix = config.get_data_prefix();
    ///
    /// assert_eq!(data_prefix, "data");
    /// ```
    pub fn get_data_prefix(&self) -> String {
        self.global
            .get("data_prefix")
            .expect("ERROR: data_prefix not found!")
            .to_string()
    }

    /// Get package name from the Config.
    ///
    /// # Returns
    ///
    /// A String containing the package name with the version attached (if any)
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let config = Config::new();
    /// let package = config.get_package_from_step(&PipelineStep::Ccs);
    ///
    /// assert_eq!(package, "pbccs/v0.1.0");
    /// ```
    pub fn get_package_from_step(&self, step: &PipelineStep) -> String {
        match step {
            PipelineStep::Minimap => {
                return step.to_str();
            }
            _ => {
                let mut package = step.to_str();

                if package == "ccs" {
                    package = String::from("pbccs");
                }

                let version = self
                    .packages
                    .get(&package)
                    .expect(format!("ERROR: Package not found -> {}", package).as_str())
                    .to_string();

                format!("{}/{}", package, version)
            }
        }
    }

    /// Get the directories for a given pipeline step.
    ///
    /// # Arguments
    ///
    /// * `step` - The pipeline step for which to get the directories.
    /// * `global_output_dir` - The global output directory.
    ///
    /// # Returns
    ///
    /// A tuple containing the input directory and the output directory for the given pipeline step.
    ///
    /// # Example
    ///
    /// ```
    /// let config = Config::new();
    /// let step = PipelineStep::Ccs;
    /// let global_output_dir = PathBuf::from("/path/to/output");
    ///
    /// let (input_dir, output_dir) = config.get_step_dirs(&step, &global_output_dir);
    ///
    /// assert_eq!(input_dir, PathBuf::from("/path/to/input"));
    /// assert_eq!(output_dir, PathBuf::from("/path/to/output/ccs"));
    /// ```
    pub fn get_step_dirs(
        &self,
        step: &PipelineStep,
        global_output_dir: &Path,
    ) -> (PathBuf, PathBuf) {
        let input_dir = self
            .get_param(*step, INPUT_DIR)
            .expect("ERROR: input_dir not found for ccs in config.toml!")
            .to_path_buf();
        let step_output_dir = global_output_dir.join(
            self.get_param(*step, OUTPUT_DIR)
                .expect("ERROR: output_dir not found for ccs in config.toml!")
                .to_path_buf(),
        );
        std::fs::create_dir_all(&step_output_dir).expect("ERROR: failed to create output dir!");

        (input_dir, step_output_dir)
    }

    /// Get custom fields for a given step.
    ///
    /// # Arguments
    ///
    /// * `step` - The step for which to retrieve custom fields.
    /// * `fields` - A vector of field names to retrieve.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let step = PipelineStep::Ccs;
    /// let config = Config::default();
    /// let custom_fields = config.get_step_custom_fields(&step, vec!["field1", "field2"]);
    ///
    /// assert_eq!(custom_fields, vec!["value1", "value2"]);
    /// ```
    pub fn get_step_custom_fields(&self, step: &PipelineStep, fields: Vec<&str>) -> Vec<String> {
        fields
            .into_iter()
            .map(|field| {
                self.get_param(*step, field)
                    .expect(
                        format!("ERROR: {} not found for {} in config.toml!", field, step).as_str(),
                    )
                    .to_string()
            })
            .collect()
    }

    /// Get arguments for a given step.
    ///
    /// # Arguments
    ///
    /// * `step` - The step for which to retrieve arguments.
    /// * `exclude` - A vector of argument names to exclude.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let step = PipelineStep::Ccs;
    /// let config = Config::default();
    /// let args = config.get_step_args(&step, vec!["arg1", "arg2"]);
    ///
    /// assert_eq!(args, "arg3 arg4");
    /// ```
    pub fn get_step_args(&self, step: &PipelineStep, exclude: Vec<&str>) -> String {
        let args = self
            .params()
            .get(step)
            .expect("ERROR: ccs not found in config.toml!")
            .flat(Some(exclude));

        args
    }

    /// Generates a unique random run ID of 4 characters.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let config = Config::default();
    /// config.set_run_id();
    /// let run_id = config.get_run_id();
    /// ```
    pub fn set_run_id(&mut self) {
        let handle = self
            .metadata
            .get_mut(RUN_ID)
            .expect("ERROR: RUN_ID not found in metadata!");

        let now = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("ERROR: Time went backwards")
            .as_nanos();

        let mut id = String::with_capacity(RUN_ID_LEN);

        // Use simple deterministic mixing to extract characters
        let mut hash = now;
        for _ in 0..RUN_ID_LEN {
            let idx = (hash % (CHARSET.len() as u128)) as usize;
            id.push(CHARSET[idx] as char);
            hash /= 7; // Crude entropy mixing
        }

        *handle = id;
    }

    /// Get the run ID.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let config = Config::default();
    ///
    /// config.set_run_id();
    /// let run_id = config.get_run_id();
    ///
    /// println!("Run ID: {}", run_id);
    /// ```
    pub fn get_run_id(&self) -> String {
        self.metadata
            .get(RUN_ID)
            .expect("ERROR: RUN_ID not found in metadata!")
            .clone()
    }

    /// Get the format package/version for a given package.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let config = Config::default();
    ///
    /// let package = "example_package";
    /// let version = config.get_package_version(package);
    ///
    /// println!("Package: {}, Version: {}", package, version);
    /// ```
    pub fn get_custom_package(&self, package: &str) -> String {
        format!(
            "{}/{}",
            package,
            self.metadata
                .get(package)
                .expect(&format!("ERROR: {} not found in config.packages!", package))
        )
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
#[derive(Deserialize, Debug, Clone, Copy, PartialEq, PartialOrd, Eq, Ord, Hash)]
pub enum PipelineStep {
    Ccs,
    Lima,
    Refine,
    Cluster,
    Minimap,
    Polya,
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
            "minimap2" => Ok(Self::Minimap),
            "polya" => Ok(Self::Polya),
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
            1 => Ok(Self::Ccs),
            2 => Ok(Self::Lima),
            3 => Ok(Self::Refine),
            4 => Ok(Self::Cluster),
            5 => Ok(Self::Minimap),
            6 => Ok(Self::Polya),
            7 => Ok(Self::LoadGenome),
            _ => Err(format!("ERROR: Invalid pipeline step: {}", i)),
        }
    }

    /// Convert a PipelineStep enum to a string.
    ///
    /// # Returns
    ///
    /// A string containing the pipeline step.
    ///
    /// # Note
    ///
    /// Both `PipelineStep::Cluster` and `PipelineStep::Minimap`
    /// are converted to their parent package 'isoseq'.
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
            Self::Refine => "isoseq".into(),
            Self::Cluster => "isoseq".into(),
            Self::Minimap => "minimap2".into(),
            Self::Polya => "polya".into(),
            Self::LoadGenome => "load-genome".into(),
        }
    }

    /// Convert a PipelineStep enum to a unique string.
    ///
    /// # Returns
    ///
    /// A string containing the pipeline step.
    ///
    /// # Note
    ///
    /// Solves `PipelineStep::Cluster` and `PipelineStep::Minimap`
    /// having 'isoseq' as their parent package.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let step = PipelineStep::Ccs;
    /// let s = step.to_str();
    ///
    /// assert_eq!(s, "ccs");
    /// ```
    pub fn to_unique_str(&self) -> String {
        match self {
            Self::Ccs => "ccs".into(),
            Self::Lima => "lima".into(),
            Self::Refine => "refine".into(),
            Self::Cluster => "cluster".into(),
            Self::Minimap => "minimap2".into(),
            Self::Polya => "polya".into(),
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
            Self::Ccs => 1,
            Self::Lima => 2,
            Self::Refine => 3,
            Self::Cluster => 4,
            Self::Minimap => 5,
            Self::Polya => 6,
            Self::LoadGenome => 7,
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

impl std::fmt::Display for PipelineStep {
    /// Format the PipelineStep enum as a string.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let step = PipelineStep::Ccs;
    ///
    /// assert_eq!(format!("{}", step), "ccs");
    /// ```
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_unique_str())
    }
}

/// A struct representing step parameters.
///
/// # Fields
///
/// * `args` - A HashMap containing step arguments.
#[derive(Deserialize, Debug, Clone)]
pub struct StepParams {
    #[serde(flatten)]
    values: HashMap<String, ParamValue>,
}

impl StepParams {
    /// Flatten the parameters into a single string for CLI execution.
    ///
    /// # Returns
    ///
    /// A string containing the flattened parameters.
    ///
    /// # Note
    ///
    /// All parameters with 2 or less characters are interpreted as short
    /// flags and will be prefixed with a single dash (`-`). Parameters with
    /// more than 2 characters will be prefixed with two dashes (`--`).
    /// Furthremore, all parameters in SPECIAL_PARAMETER will be suffixed with an
    /// equal sign (`=`).
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let params = StepParams {
    ///    values: HashMap::new(),
    /// };
    ///
    /// let flat = params.flat();
    ///
    /// assert_eq!(flat, "");
    /// ```
    pub fn flat(&self, exclude: Option<Vec<&str>>) -> String {
        let exclude = exclude
            .unwrap_or_default()
            .into_iter()
            .collect::<HashSet<_>>();

        self.values
            .iter()
            .filter(|(key, _)| !exclude.contains(key.as_str()))
            .map(|(key, value)| {
                let mut argument = if key.len() > 2 {
                    if SPECIAL_PARAMETER.contains(&key.as_str()) {
                        format!("--{}=", key)
                    } else {
                        format!("--{} ", key)
                    }
                } else {
                    if SPECIAL_PARAMETER.contains(&key.as_str()) {
                        format!("-{}=", key)
                    } else {
                        format!("-{} ", key)
                    }
                };

                match value {
                    ParamValue::Int(i) => argument.push_str(&i.to_string()),
                    ParamValue::Float(flt) => argument.push_str(&flt.to_string()),
                    ParamValue::Bool(b) => argument.push_str(&b.to_string()),
                    ParamValue::Str(s) => argument.push_str(&s.clone()),
                }

                argument
            })
            .collect::<Vec<_>>()
            .join(" ")
    }

    /// Get a parameter value from a StepParams struct.
    ///
    /// # Arguments
    ///
    /// * `key` - A string containing the parameter key.
    ///
    /// # Returns
    ///
    /// An Option containing the parameter value.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let params = StepParams {
    ///   values: HashMap::new(),
    /// };
    ///
    /// let value = params.get("min-rq");
    ///
    /// assert_eq!(value, None);
    /// ```
    pub fn get(&self, key: &str) -> Option<&ParamValue> {
        self.values.get(key)
    }
}

/// Represents a parameter value for any step
///
/// # Example
///
/// ``` rust, no_run
/// let value = ParamValue::Int(1);
///
/// assert_eq!(value, ParamValue::Int(1));
/// ```
#[derive(Deserialize, Debug, Clone)]
#[serde(untagged)]
pub enum ParamValue {
    Int(i64),
    Float(f64),
    Bool(bool),
    Str(String),
}

impl ParamValue {
    /// Convert a ParamValue to a PathBuf.
    ///
    /// # Returns
    ///
    /// A PathBuf.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let value = ParamValue::Str("path/to/file".into());
    ///
    /// assert_eq!(value.to_path_buf(), PathBuf::from("path/to/file"));
    /// ```
    pub fn to_path_buf(&self) -> PathBuf {
        match self {
            ParamValue::Str(s) => PathBuf::from(s),
            _ => PathBuf::new(),
        }
    }

    /// Convert a ParamValue to a string.
    ///
    /// # Returns
    ///
    /// A string.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let value = ParamValue::Str("string".into());
    ///
    /// assert_eq!(value.to_string(), "string");
    /// ```
    pub fn to_string(&self) -> String {
        match self {
            ParamValue::Str(s) => s.clone(),
            ParamValue::Int(i) => i.to_string(),
            ParamValue::Float(f) => f.to_string(),
            ParamValue::Bool(b) => b.to_string(),
        }
    }

    /// Convert a ParamValue to an integer.
    ///
    /// # Returns
    ///
    /// An integer.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let value = ParamValue::Int(1);
    ///
    /// assert_eq!(value.to_int(), 1);
    /// ```
    pub fn to_int(&self) -> i64 {
        match self {
            ParamValue::Int(i) => *i,
            _ => 0,
        }
    }

    /// Convert a ParamValue to a float.
    ///
    /// # Returns
    ///
    /// A float.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let value = ParamValue::Float(1.0);
    ///
    /// assert_eq!(value.to_float(), 1.0);
    /// ```
    pub fn to_float(&self) -> f64 {
        match self {
            ParamValue::Float(f) => *f,
            _ => 0.0,
        }
    }

    /// Convert a ParamValue to a boolean.
    ///
    /// # Returns
    ///
    /// A boolean.
    ///
    /// # Example
    ///
    /// ``` rust, no_run
    /// let value = ParamValue::Bool(true);
    ///
    /// assert_eq!(value.to_bool(), true);
    /// ```
    pub fn to_bool(&self) -> bool {
        match self {
            ParamValue::Bool(b) => *b,
            _ => false,
        }
    }
}

/// Implement Display for ParamValue
///
/// # Example
///
/// ``` rust, no_run
/// let value = ParamValue::Int(1);
///
/// assert_eq!(format!("{}", value), "1");
/// ```
impl std::fmt::Display for ParamValue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParamValue::Int(i) => write!(f, "{}", i),
            ParamValue::Float(flt) => write!(f, "{}", flt),
            ParamValue::Bool(b) => write!(f, "{}", b),
            ParamValue::Str(s) => write!(f, "{}", s),
        }
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

/// Deserialize a HashMap of PipelineStep enums and StepParams.
///
/// # Arguments
///
/// * `deserializer` - A serde Deserializer.
///
/// # Returns
///
/// A Result containing a HashMap of PipelineStep enums and StepParams or an error.
///
fn deserialize_to_hash<'de, D>(
    deserializer: D,
) -> Result<HashMap<PipelineStep, StepParams>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    let raw: HashMap<String, StepParams> = HashMap::deserialize(deserializer)?;

    raw.into_iter()
        .map(|(key, value)| PipelineStep::from_str(&key).map(|step| (step, value)))
        .collect::<Result<HashMap<_, _>, _>>()
        .map_err(serde::de::Error::custom)
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

/// Executes a shell command and logs the output.
///
/// # Arguments
///
/// * `cmd` - The shell command to execute.
/// * `log_msg` - The message to log upon successful execution.
/// * `tool` - The name of the tool being executed.
///
/// # Example
///
/// ```rust, no_run
/// shell("ls -l".to_string(), "Listing files", "");
/// ```
pub fn shell(cmd: String, log_msg: &str, tool: &str) {
    let tool = if tool.is_empty() { ISOPIPE } else { tool };

    let output = std::process::Command::new("sh")
        .arg("-c")
        .arg(cmd.clone())
        .output()
        .expect("ERROR: Failed to execute process");

    if output.status.success() {
        log::info!("INFO [{}]: {}!", tool, log_msg);
    } else {
        log::error!(
            "ERROR: failed to execute {}\n{}",
            cmd,
            String::from_utf8_lossy(&output.stderr)
        );
        std::process::exit(1);
    }
}
