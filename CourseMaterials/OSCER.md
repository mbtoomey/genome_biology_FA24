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

![WinScp](https://github.com/mbtoomey/genome_biology_FA24/blob/main/CourseMaterials/winscp.png)

## Your directories on OSCER

* **Home** `/home/USER_ID/`
    * Persistent storage - will last as long as you have an active account
    * limited to 20 GB
    * I mostly use this for scripts and smaller template files like assembled genomes or transcriptomes that I am aligning raw reads to. 
    
* **Scratch** `/scratch/USER_ID?`
    * Large-scale :warning: **Temporary** space
    * Files are deleted ~every two weeks
    

    




