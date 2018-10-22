This step-by-step instruction presumes that you are using AWS with FPGA Developer AMI. We strongly recommend you using this AMI. But if you insist to work on a different Linux distribution or platform, please refer to its manual to install Git, essential building tools and OpenMPI, and to compile the Graph500 code.

Please be noted that we have a special version of Graph500 without BFS kernel, and with a modified SSSP kernel provided on https://github.com/UCLA-VAST/cs-259-graph500-sssp. You can start with this.

You do need to install OpenMPI library to build the code. But you don’t need to understand it since we have removed all the MPI related code from the kernel. You can just optimize this part as if it is a normal single-process application. Depend on your distribution, you may use a different command.

To install and load OpenMPI library on your system:

    sudo yum -y install openmpi-devel
    module load /etc/modulefiles/mpi/openmpi-x86_64

Please make sure you load the modulefile provided by OpenMPI each time you log in. The path of the modulefile `/etc/modulefiles/mpi/openmpi-x86_64` may vary on a different distribution. You can check the manual or `/etc/modulefiles/` folder for the exact path.

You can checkout the code by:

    git clone https://github.com/UCLA-VAST/cs-259-graph500-sssp.git

If you see an error message, please make sure you have Git installed on your system.

To build and run the code:

    module load /etc/modulefiles/mpi/openmpi-x86_64 # required on each login
    cd graph500_project/src/
    make
    ./graph500_sssp 18 # change 18 to your SCALE

## Acceleration Evaluation

Your task is to optimize the single-source shortest path kernel in `src/sssp_reference.c`. Please submit the result of the following two commands along with your source code.

    ./graph500_sssp 18
    ./graph500_sssp 23

You should refer to the line with “sssp harmonic_mean_TEPS” as your final performance metrics, and make sure your code passes all the validation. For example:

    …
    Running SSSP 63
    Time for SSSP 63 is 0.184641
    TEPS for SSSP 63 is 2.27115e+07
    Validating SSSP 63
    Validate time for SSSP 63 is 0.089011
    …
    sssp min_TEPS:                  1.6315e+07
    sssp firstquartile_TEPS:        2.05992e+07
    sssp median_TEPS:               2.25539e+07
    sssp thirdquartile_TEPS:        2.42411e+07
    sssp max_TEPS:                  3.43145e+07
    sssp harmonic_mean_TEPS:     !  2.23902e+07
    sssp harmonic_stddev_TEPS:      411217
    sssp min_validate:              0.0851015
    sssp firstquartile_validate:    0.0858768
    sssp median_validate:           0.0869632
    sssp thirdquartile_validate:    0.0880603
    sssp max_validate:              0.0893905
    sssp mean_validate:             0.0870664
    sssp stddev_validate:           0.00129654

Please do not modify the validation code.
