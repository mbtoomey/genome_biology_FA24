# Working with the OSCER

## Logging In

There are multiple "login nodes" schooner1-3. One of these will be automatically assigned when you login to `schooner.oscer.ou.edu`. To login simple type: 

```
ssh [YOUR_ACCOUNT_ID]@schooner.oscer.ou.edu
```
into you terminal application. This will then prompt you for your password. *note your password will be hidden from the screen and nothing will appear at the cursor when you enter your password

There is also a specific node to use for large data transfers - `dtn2.oscer.ou.edu`. You should use this if you are uploading or downloading large amounts of data (~>2 GB).

## Login nodes vs. Compute Nodes

**Login node**

* Primary uses - launching and monitoring experiments via SLURM (see below)
* Can be used for limited file manipulation and small debugging steps
* :skull_and_crossbones: **DO NOT RUN ANALYSES ON THIS NODE**. Computationally intensive jobs will break the node and disrupt all users on that node. 

**Compute node**

* This is where the real work is done
* We will access these by submitting jobs (via .sbatch files) to the SLURM scheduler.

## Uploading and downloading from OSCER

I find the easiest way to move files to and from my PC to the OSCER system is with [WinScp](https://winscp.net/eng/downloads.php) transfer client. 

To link to OSCER use SFTP protocol, port 22, host: `schooner.oscer.ou.edu` or `dtn2.oscer.ou.edu` for large transfers, and enter your ID and pasword. For example:  

![WinScp](https://github.com/mbtoomey/genome_biology_FA24/blob/main/CourseMaterials/winscp.png =250x)

## Your directories on OSCER

* **Home** `/home/USER_ID/`
    * Persistent storage - will last as long as you have an active account
    * limited to 20 GB
    * I mostly use this for scripts and smaller template files like assembled genomes or transcriptomes that I am aligning raw reads to. 
    
* **Scratch** `/scratch/USER_ID?`
    * Large-scale :warning: **Temporary** space
    * Files are deleted ~every two weeks
    
## Submitting jobs to the compute node

Jobs on the OSCER are managed through [SLURM](https://slurm.schedmd.com/overview.html) a system that tracks computing resources and allocates them among requests. 

To run a job: 

* submit a job to the appropriate partition with a .sbatch file (details below)
* When a compute node becomes available and you are at the front of the queue, the node will pick up the job and begin execution
* While the job is executed output and error meaasages will be written to log and error files you specify in the .sbatch file
* Once the job is complete it is removed from the queue. 

### Job Partitions/Queues on Schooner

You will submit your job to a specific partition depending on the resources you need and the timeline of your job. As a general rule, the more resource and time you request the loger you will have to wait for your job to run. 

A few of the partitions available: 

* debug_5min: quick pick up; cannot run for more than 5 minutes
* debug: relatively quick pick up; 30 min limit
* normal: where many experiments will be; 48 hour limit
   * up to 30 GB of RAM per node
   * up to 20 cores per node
* large_mem: used when memory needs exceed the normal nodes; 48 hour limit. You will not need this for any of the course activities, but might need this for intensive jobs in your own research (e.g. de novo transcriptome assembly)
   * up to 1000 GB of RAM per node
   * up to 30 cores per nodes
   
### Our sbatch template

```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem 16G
#SBATCH --output=[job_name]_%J_stdout.txt
#SBATCH --error=[job_name]_%J_stderr.txt
#SBATCH --job-name=[job_name]
# 

bash [path to your .sh script]
```

* `#SBATCH --partition=normal` directs the job to the partition of your choice
* `#SBATCH --ntasks=1` We will always set this to "1", however someday in the future you might have an application allows for parallel process, more details [here](https://www.ou.edu/oscer/support/running_jobs_schooner)
* `#SBATCH --cpus-per-task=1` If your application allows you to use multiple processors/cores to speed up processing you will need to request them here. 
* `#SBATCH --mem 16G` Here you will set the amount of RAM memory needed for your job. Often the software documentation will give recommendations for miniumums
* `#SBATCH --output=[job_name]_%J_stdout.txt` this specifies the name of the output file. This is where any of the information that would typically appear in the console, will be written. These files will be save in the directory that you are in when you submit your job. 
* `#SBATCH --error=[job_name]_%J_stderr.txt` this specifies the name of the error log file. This is where any error messages will be written. 
* `#SBATCH --job-name=[job_name]` this is where you specify the name of your job that will appear in the SLURM queue. 

You can download the sbatch template here: [template.sbatch](https://github.com/mbtoomey/genome_biology_FA24/blob/main/CourseMaterials/template.sbatch)



    




