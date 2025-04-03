use crate::config::PipelineStep;

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
            PipelineStep::Refine => "isoseq3 refine",
            PipelineStep::Cluster => "isoseq3 cluster",
            PipelineStep::FilterQuality => "isotools iso-polya filter",
            PipelineStep::Minimap => "minimap3",
            PipelineStep::LoadGenome => "load_genome",
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
    pub fn arg(mut self, arg: &str) -> Self {
        self.cmd.push(' ');
        self.cmd.push_str(arg);
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
