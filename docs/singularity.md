# CarrierSeq Singularity
Singularity is a container, similar to Docker, that is secure to run in HPC environments. For this pipeline, we provide two example containers, each of which provides modular access to the different software inside. By way of using the Standard Container Integration Format (SCI-F) with Singularity, we have a lot of freedom in deciding on what level of functions we want to expose to the user. A developer will want easy access to the core tools (e.g., bwa, seqtk) while a user will likely want one level up, on the level of a collection of steps associated with some task (e.g., mapping). We will walk through the steps of building and using each one.

## Setup
You first need to install Singularity. For this tutorial, we are using the development branch with the upcoming version 2.4. You can install it doing the following.

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
sudo singularity build --writable carrierseq.img Singularity
```

Note that `--writable` is only being done because we potentially want to develop, or
tweak our image. Once we are happy and want to build "for reals" we would run this:

```
sudo singularity build carrierseq.img Singularity
```

This will build an essentially frozen image, so your pipeline is forever preserved.

## Listing Pipeline Apps
If you didn't know anything about the image, you would want to explore it. SCI-F 
provides easy commands to do that. 


### Global Help
```
singularity help carrierseq.img

    CarrierSeq is a sequence analysis workflow for low-input nanopore 
    sequencing which employs a genomic carrier.

    Github Contributors: Angel Mojarro (@amojarro), 
                         Srinivasa Aditya Bhattaru (@sbhattaru),   
                         Christopher E. Carr (@CarrCE),
                         Vanessa Sochat (@vsoch).

    fastq-filter from: https://github.com/nanoporetech/fastq-filter
    see:
           singularity run --app readme carrierseq.img | less for more detail
   
```

The command mentioned at the body is an internal provided only to add, preserve, and
make the README accessible. Cool!

```
singularity run --app readme carrierseq.img | less

# CarrierSeq

## About

bioRxiv doi: https://doi.org/10.1101/175281

CarrierSeq is a sequence analysis workflow for low-input nanopore sequencing which employs a genomic carrier.

Github Contributors: Angel Mojarro (@amojarro), Srinivasa Aditya Bhattaru (@sbhattaru), Christopher E. Carr (@CarrCE), and Vanessa Sochat (@vsoch).</br> 
fastq-filter from: https://github.com/nanoporetech/fastq-filter

.......etc
```

We can see all apps in the image:

```
singularity apps carrierseq.img
mapping
poisson
readme
sorting
sra-tooklit
```

And then we can ask for help for any of the pipeline steps:

```
singularity help --app mapping carrierseq.img
singularity help --app poisson carrierseq.img
singularity help --app sorting carrierseq.img
```

We can also look at metadata for the entire image, or for an app.  The inspect
command can expose environment variables, labels, the definition file, tests, and
runscripts.

See `singularity inspect --help` for more examples:

```
singularity inspect carrierseq.img

{
    "org.label-schema.usage.singularity.deffile.bootstrap": "docker",
    "org.label-schema.usage.singularity.deffile": "Singularity",
    "org.label-schema.usage": "/.singularity.d/runscript.help",
    "org.label-schema.schema-version": "1.0",
    "org.label-schema.usage.singularity.deffile.from": "ubuntu:14.04",
    "org.label-schema.build-date": "2017-09-20T18:16:50-07:00",
    "BIORXIV_DOI": "https://doi.org/10.1101/175281",
    "org.label-schema.usage.singularity.runscript.help": "/.singularity.d/runscript.help",
    "org.label-schema.usage.singularity.version": "2.3.9-development.gaaab272",
    "org.label-schema.build-size": "1419MB"
}


singularity inspect --app mapping carrierseq.img
{
    "FQTRIM_VERSION": "v0.9.5",
    "SEQTK_VERSION": "v1.2",
    "BWA_VERSION": "v0.7.15",
    "SINGULARITY_APP_NAME": "mapping",
    "SINGULARITY_APP_SIZE": "9MB"
}
```

## Overall Strategy
Since the data is rather large, we are going to map a folder to our $PWD where the analysis is to run. This directory, just like the modular applications, has a known and predictable location. So our steps are going to look like this:

```
# 0. Make an empty random folder to bind to for data.
mkdir data

# 1. Download data, map the data base to an empty folder on our local machine
singularity run --app sra-toolkit --bind data:/scif/data carrierseq.img

# 2. Perform mapping step of pipeline, mapping the same folder.
singularity run --app mapping --bind data:/scif/data carrierseq.img

# 3. perform poisson regression on filtered reads
singularity run --app poisson --bind data:/scif/data carrierseq.img

# 4. Finally, sort results
singularity run --app sorting --bind data:/scif/data carrierseq.img
```

For any of the above steps, see the `singularity help --app <appname> carrierseq.img` for how to
customize settings and environment variables. This demo is intended to use the defaults.


## CarrierSeq Pipeline
The common user might want access to the components of the pipeline that the software provides. In the case of CarrierSeq, this means mapping, poisson, and then sorting. If the image isn't provided (e.g., a Singularity Registry or Singularity Hub) the user can build from the build recipe file, `Singularity`:

```
sudo singularity build carrierseq.img Singularity
```

As a user, you want a container that exposes enough metadata to run different steps of the pipeline, but you don't want to need to know the specifics of each command call or path within the container. 
## CarrierSeq Development


For the 
1. Download and install Singularity (https://singularityware.github.io/install-linux). For this demo, we are using the development branch (2.4)
2. Build `sudo singularity create --size 4000 carrierseq.img && sudo singularity bootstrap carrierseq.img Singularity`
3. Run container!

