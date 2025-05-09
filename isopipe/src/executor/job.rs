use crate::config::PipelineStep;
use std::fmt::Write;

/// Struct to represent a job to be executed
/// by the pipeline
///
/// # Example
///
/// ```rust, no_run
/// use isopipe::executor::job::Job;
///
/// let job = Job::new()
///    .task(PipelineStep::Ccs)
///     .arg("input.bam")
///     .arg("output.bam")
///     .arg("chunks");
///
/// assert_eq!(job.cmd, "ccs input.bam output.bam chunks");
/// ```
#[derive(Debug, Clone)]
pub struct Job {
    pub cmd: String,
}

impl Job {
    /// Create a new job
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use isopipe::executor::job::Job;
    ///
    /// let job = Job::new();
    ///
    /// assert_eq!(job.cmd, "");
    /// ```
    pub fn new() -> Self {
        Self { cmd: String::new() }
    }

    /// Create a new job from a command string
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use isopipe::executor::job::Job;
    ///
    /// let job = Job::from("ccs input.bam output.bam chunks");
    ///
    /// assert_eq!(job.cmd, "ccs input.bam output.bam chunks");
    /// ```
    pub fn from(cmd: String) -> Self {
        Self { cmd }
    }

    /// Add a task to the job
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use isopipe::executor::job::{Job, PipelineStep};
    ///
    /// let job = Job::new()
    ///    .task(PipelineStep::Ccs);
    ///
    /// assert_eq!(job.cmd, "ccs");
    /// ```
    pub fn task(mut self, step: PipelineStep) -> Self {
        let step_cmd = match step {
            PipelineStep::Ccs => "ccs",
            PipelineStep::Lima => "lima",
            PipelineStep::Refine => "isoseq refine",
            PipelineStep::Cluster => "isoseq cluster",
            PipelineStep::Polya => "",
            PipelineStep::Minimap => "minimap2",
            PipelineStep::Fusion => "isotools iso-fusion",
            PipelineStep::Orf => "",
        };

        self.cmd.push_str(step_cmd);
        self
    }

    /// Add an argument to the job
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use isopipe::executor::job::Job;
    ///
    /// let job = Job::new()
    ///     .task(PipelineStep::Ccs)
    ///     .arg("input.bam");
    ///
    /// assert_eq!(job.cmd, "ccs input.bam");
    /// ```
    pub fn arg<D: std::fmt::Display>(mut self, arg: D) -> Self {
        self.cmd.push(' ');
        write!(&mut self.cmd, "{arg}").expect("ERROR: Failed to append arg to cmd!");

        self
    }

    /// Add multiple arguments to the job
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use isopipe::executor::job::Job;
    ///
    /// let job = Job::new()
    ///     .task(PipelineStep::Ccs)
    ///     .args(&["input.bam", "output.bam", "chunks"]);
    ///
    /// assert_eq!(job.cmd, "ccs input.bam output.bam chunks");
    /// ```
    pub fn args(mut self, args: &[&str]) -> Self {
        for arg in args {
            self.cmd.push(' ');
            self.cmd.push_str(arg);
        }
        self
    }

    pub fn cmd(&self) -> &str {
        &self.cmd
    }
}
