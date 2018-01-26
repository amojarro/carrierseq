# CarrierSeq Singularity
Singularity is a container, similar to Docker, that is secure to run in HPC environments. For this pipeline, we provide two example containers, each of which provides modular access to the different software inside. By way of using a Scientific Filesystem (SCIF) with Singularity, we have a lot of freedom in deciding on what level of functions we want to expose to the user. A developer will want easy access to the core tools (e.g., bwa, seqtk) while a user will likely want one level up, on the level of a collection of steps associated with some task (e.g., mapping). We will walk through the steps of building and using each one.

## Setup
You first need to install Singularity. For this tutorial, we are using the development branch with the current version 2.4.x. You can install it doing the following.

```
git clone -b development https://www.github.com/singularityware/singularity.git
cd singularity
./autogen.sh
./configure --prefix=/usr/local
make
sudo make install
```

Now you are ready to go!

## Get Data
If you aren't familar with genomic analysis (as I'm not) you will be utterly confused about
how to download data. The reference is [provided](../reference), but the input data isn't. The good news is, the tool to download the data is part of the container. So you don't need to do the nonsense that I did to get it. 

## Build the image
Building looks like this:

```
sudo singularity build cseq Singularity
```

This will build a read only (essentially frozen) image, so your pipeline is forever preserved.


## Exploring the Container
If you didn't know anything about the image, you would want to explore it. SCIF 
provides easy commands to do that. When we just execute the container, the SCIF shows us
it's default help:

```
./cseq 
[help] executing /bin/bash /scif/apps/help/scif/runscript
    CarrierSeq is a sequence analysis workflow for low-input nanopore
    sequencing which employs a genomic carrier.
    Github Contributors: Angel Mojarro (@amojarro),
                         Srinivasa Aditya Bhattaru (@sbhattaru),
                         Christopher E. Carr (@CarrCE),
                         Vanessa Sochat (@vsoch).
    fastq-filter from: https://github.com/nanoporetech/fastq-filter
    To see applications installed in the Scientific Filesystem:
    scif apps
    To run a typical pipeline, you might do:
    scif run mapping
    scif run poisson
    scif run sorting
    If you install in a container, the entrypoint should be scif, and then
    issue the above commands to it.
```

### SCI-F Apps

You can see apps in the image, as instructed:

```
./cseq apps
  download
      help
   mapping
   poisson
    readme
 reference
   sorting
```

or see help for a specific app. For example, readme exists only to capture and then
print the entire readme of the repository (or any repository container README.md for
which the SCIF is installed for) to the console:

```
./cseq run readme | tail
The matrix illustrates the reads/channel distribution of B. subtilis, contamination, and HQNRs across all 512 nanopore channels. Here we are able to visually identify overly productive channels (e.g., 191 reads/channel, etc) producing likely HQNRs.
![alt text](https://github.com/amojarro/carrierseq/blob/master/example/carrierseq_roi_q9_p005.png)

### HQNR Pore Occupancy
“Bad” channels identified by CarrierSeq as HQNR-associated (reads/channel > 7).
![alt text](https://github.com/amojarro/carrierseq/blob/master/example/carrierseq_hqnrs_q9_p005.png)

### Target Reads Pore Occupancy
“Good” channels identified by CarrierSeq as non-HQNR-associated (reads/channel ≤ 7). Channels producing 6 or more reads yield HQNRs that have satisfied our CarrierSeq parameters. By imposing a stricter p-value, CarrierSeq may be able to reject more HQNRs (e.g., xcrit = 5).
![alt text](https://github.com/amojarro/carrierseq/blob/master/example/carrierseq_target_reads_q9_p005.png)
```

And then we can ask for help for any of the pipeline steps:

```
./cseq help download
```

We can also look at metadata for the entire image, or for an app.  The inspect
command can expose environment variables, labels, the definition file, tests, and
runscripts.

```
# Returns a large data structure with all apps
./cseq inspect

# Returns for a particular app
./cseq inspect download
{
    "download": {
        "apphelp": [
            "   The original sra-toolkit does not serve the correct data, so for now you ",
            "   should download data from ",
            "   https://www.dropbox.com/sh/vyor82ulzh7n9ke/AAC4W8rMe4z5hdb7j4QhF_IYa?dl=0) and then move into some data folder you intend to mount:",
            "      mv $HOME/Downloads/all_reads.fastq data/"
        ]
    }
}
```

## Overall Strategy
Since the data is rather large, we are going to map a folder to our $PWD where the analysis is to run. This directory, just like the modular applications, has a known and predictable location. So our steps are going to look like this:

```
# 0. Make an empty random folder to bind to for data.
mkdir data

# 1. Download data, map the data base to an empty folder on our local machine
#     (we actually will do this from the browser as the sri toolkit is changed.
# 2. Perform mapping step of pipeline, mapping the same folder.
# 3. perform poisson regression on filtered reads
# 4. Finally, sort results

# Download data from https://www.dropbox.com/sh/vyor82ulzh7n9ke/AAC4W8rMe4z5hdb7j4QhF_IYa?dl=0
# See issue --> https://github.com/amojarro/carrierseq/issues/1
cseq="singularity run --bind data:/scif/data cseq"

# $cseq run <app>
```

## CarrierSeq Pipeline
The common user might want access to the larger scoped pipeline that the software provides. In the case of CarrierSeq, this means (optionally, download) mapping, poisson, and then sorting. If the image isn't provided (e.g., a Singularity Registry or Singularity Hub) the user can build from the build recipe file, `Singularity.devel`. If you haven't done this already:

```
sudo singularity build cseq Singularity.devel
```

The author is still working on updating the automated download, for now download from [here](https://www.dropbox.com/sh/vyor82ulzh7n9ke/AAC4W8rMe4z5hdb7j4QhF_IYa?dl=0) and then move into some data folder you intend to mount:

``
mv $HOME/Downloads/all_reads.fastq data/
```

Let's be lazy and put the bind and singularity command into an environment variable so we don't need
to type it many times.

```
cseq="singularity run --bind data:/scif/data cseq"
```

And then the various steps of the pipeline are run as was specified above:

```
$cseq run mapping
$cseq run poisson
$cseq run sorting
```

There are even helper functions that the main ones use to get various environments. For example, the `reference` app serves only to return the path to the reference genome!

```
$ ./cseq exec reference echo [e]reference_genome
[reference] executing /bin/echo $reference_genome
$SCIF_APPROOT_reference/lambda_ecoli.fa
```

or run quietly

```
$ ./cseq --quiet exec reference echo [e]reference_genome
$SCIF_APPROOT_reference/lambda_ecoli.fa
```
and what is that evaluated to?

```
./cseq --quiet exec reference echo [e]SCIF_APPROOT_reference/lambda_ecoli.fa
/scif/apps/reference/lambda_ecoli.fa
```

Cool!


### 1. Mapping
To run mapping, bind the data folder to the image, and specify the app to be mapping:

```
singularity run --bind data:/scif/data cseq run mapping
```

## 2. Poisson

```
singularity run --bind data:/scif/data cseq run poisson
```

## 3. Sorting

```
singularity run --bind data:/scif/data cseq run sorting
```


## How Can I Change It?
Given two containers that share the same input and output schema, I could
swap out of the steps:


```
...
singularity run --bind data:/scif/data container1 run sorting
singularity run --bind data:/scif/data container2 run sorting
```

or even provide a single container with multiple options for the same step


```
...
singularity run --bind data:/scif/data cseq run sorting1
singularity run --bind data:/scif/data cseq run sorting2
```

As a user, you want a container that exposes enough metadata to run different steps of the pipeline, but you don't want to need to know the specifics of each command call or path within the container. In the above, I can direct the container to my mapped input directory
and specify a step in the pipeline, and I dont need to understand how to use `bwa` or `grep` or `seqtk`, or any of the other software
that makes up each.


### Interactive Environment
If you are curious about all the Scientific Filesystem variables available to you:

```
$cseq exec mapping env | grep SCIF_
```

and read [the specification](https://sci-f.github.io/specification) for a complete list.

Or if you wanted an interactive shell to explore running commands inside the container:

```
singularity shell --bind data:/scif/data cseq
```

or enter a shell in context of a scif application (mapping)

```
$ singularity run cseq shell mapping
$ echo $SCIF_APPNAME
mapping
```

## CarrierSeq Development
The developer has a different use case - to have easy command line access to the lowest level of executables installed in the container. Given a global install of all software, without SCI-F I would need to look at `$PATH` to see what has been added to the path, and then list executables in path locations to find new software installed to, for example, `/usr/bin`. There is no way to easily and programatically "sniff" a container to understand the changes, and the container is a black development box, perhaps only understood by the creator or with careful inspection of the build recipe.

In this use case, we want to build the development container.

```
sudo singularity build carrierseq.dev.img Singularity.devel
```

Now when we look at apps, we see all the core software that can be combined in specific ways to produce a pipeline step like "mapping".

```
singularity apps carrierseq.dev.img
bwa
fqtrim
python
seqtk
sra-toolkit
```

each of which might be run, exec to activate the app environment, or shell to shell into the container under the context of a specific app:

```
# Open interactive python
singularity run --app python carrierseq.dev.img

# Load container with bwa on path
$ singularity shell --app bwa carrierseq.dev.img
$ which bwa
$ /scif/apps/bwa/bin/bwa
```

These two images that serve equivalent software is a powerful example of the flexibility of SCI-F. The container creator can choose the level of detail to expose to a user that doesn't know how it was created, or perhaps has varying levels of expertise. A lab that is using core tools for working with sequence data might have preference for the development container, while a finalized pipeline distributed with a publication would have preference for the first. In both cases, the creator doesn't need to write custom scripts for a container to run a particular app, or to expose environment variables, tests, or labels. By way of using SCI-F, this happens automatically. 
