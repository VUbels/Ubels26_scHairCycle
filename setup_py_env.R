###############################################################
# SETUP CONDA PY ENVIRONMENT WITH APPROPRIATE GPU COMPATABILITY
###############################################################


setup_py_env <- function(py_env_name, py_location) {
  
  # PYTORCH_ROCM_ARCH type necessary from compatability matrix if using AMD ROCm
  # Check https://rocm.docs.amd.com/projects/install-on-linux/en/latest/reference/system-requirements.html

  Sys.setenv(PYTORCH_ROCM_ARCH = "gfx1201")
  options(reticulate.conda_binary = py_location)
  
  envs <- reticulate::conda_list()$name
  
  if (!(py_env_name %in% envs)) {
    
    message("Creating new environment: ", py_env_name)
    
    # Create env
    
    reticulate::conda_create(
      envname = py_env_name,
      python_version = "3.12",
      packages = c("pip", "umap-learn", "tables", "scanorama")
    )
    
    # Activate
    
    reticulate::use_condaenv(py_env_name, required = TRUE)
    
    # Get pip path for this environment
    conda_envs <- reticulate::conda_list()
    env_info <- conda_envs[conda_envs$name == py_env_name, ]
    pip_path <- file.path(dirname(dirname(env_info$python)), "bin", "pip")
    
    # Install PyTorch ROCm using index URL, change for CUDA/ROCm (Nvidia/AMD)
    
    message("Installing PyTorch ROCm...")
    system(paste(pip_path, "install torch torchvision --index-url https://download.pytorch.org/whl/rocm6.4"))
    
    # Remove bundled HSA runtime specific for WSL2/Linux/ROCm to ensure GPU can be found
    # For more info see: https://github.com/ROCm/ROCm/issues/4682\
    
    message("Applying HSA runtime fix...")
    torch_lib_path <- reticulate::py_capture_output({
      reticulate::py_run_string("import torch; import os; print(os.path.join(os.path.dirname(torch.__file__), 'lib'))")
    })
    torch_lib_path <- trimws(torch_lib_path)
    system(paste("rm -f", file.path(torch_lib_path, "libhsa-runtime64.so*")))
    
    # Install CellBender with custom pull 420 to make python 3.12/numpy 2.0+ compatible
    
    message("Installing CellBender...")
    system(paste(pip_path, "install git+https://github.com/broadinstitute/CellBender.git@refs/pull/420/head"))
    
    message("\n###########################")
    message("Environment setup complete!")
    message("R session needs to restart for GPU detection to work.")
    message("#############################\n")
    
    # Auto-restart if in RStudio, otherwise prompt
    
    if (Sys.getenv("RSTUDIO") == "1") {
      message("Restarting RStudio session in 3 seconds...")
      Sys.sleep(3)
      .rs.restartR()
    } else {
      message("Please restart R manually, then re-run your script.")
      message("Type: quit() or use your IDE's restart option")
      return(NULL)
    }
    
  } else {
    # Environment exists - verify it works
    
    reticulate::use_condaenv(py_env_name, required = TRUE)
    torch <- reticulate::import("torch")
    gpu_available <- torch$cuda$is_available()
    
    conda_envs <- reticulate::conda_list()
    env_info <- conda_envs[conda_envs$name == py_env_name, ]
    python_path <- env_info$python
    cellbender_bin <- file.path(dirname(dirname(env_info$python)), "bin", "cellbender")
    
    if (gpu_available) {
      message("Using existing env with GPU: ", torch$cuda$get_device_name(0L))
    } else {
      warning("Existing env - GPU not detected!")
    }
    
    return(list(
      gpu_ok = gpu_available,
      cellbender_bin = cellbender_bin,
      python_path = python_path 
    ))
  }
}