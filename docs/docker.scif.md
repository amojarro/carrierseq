# CarrierSeq Scientific Filesystem (Docker)
The Singularity [container described](singularity.md) works by way of installing a Scientific Filesystem, meaning that each of the different functions in the pipeline is an internally modular application in the container. We can install these same applications to a Docker image in the case that we have preference for this technology. This set of documentation will walk through the steps to do that.

## Setup
Follow the instructions in the primary [README.md](../README.md) for installing Docker. Then download data. The reference is [provided](../reference), but the input data isn't. Instructions for getting this data will be provided in the tutorial below.

## Build the image
Building looks like this. Notice that we aren't using the primary [Dockerfile](../Dockerfile) in the repository, we are using the [Dockerfile.scif](../Dockerfile.scif) that is the version for the Scientific Filesystem.

```
docker build -f Dockerfile.scif -t vanessa/cseq .
```

## Exploring the Container
If you didn't know anything about the image, you would want to explore it. SCIF 
provides easy commands to do that. When we just execute the container, the SCIF shows us

Let's define a quicker entry point to using the container:

```
cseq="docker run vanessa/cseq"
```

### SCI-F Apps

You can see apps in the image:

```
$cseq apps
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
$cseq run readme | tail
The matrix illustrates the reads/channel distribution of B. subtilis, contamination, and HQNRs across all 512 nanopore channels. Here we are able to visually identify overly productive channels (e.g., 191 reads/channel, etc) producing likely HQNRs.
![alt text](https://github.com/amojarro/carrierseq/blob/master/example/carrierseq_roi_q9_p005.png)

### HQNR Pore Occupancy
"Bad" channels identified by CarrierSeq as HQNR-associated (reads/channel > 7).
![alt text](https://github.com/amojarro/carrierseq/blob/master/example/carrierseq_hqnrs_q9_p005.png)

### Target Reads Pore Occupancy
"Good" channels identified by CarrierSeq as non-HQNR-associated (reads/channel ≤ 7). Channels producing 6 or more reads yield HQNRs that have satisfied our CarrierSeq parameters. By imposing a stricter p-value, CarrierSeq may be able to reject more HQNRs (e.g., xcrit = 5).
![alt text](https://github.com/amojarro/carrierseq/blob/master/example/carrierseq_target_reads_q9_p005.png)
```

And then we can ask for help for any of the pipeline steps:

```
$cseq help download
```

We can also look at metadata for the entire image, or for an app.  The inspect
command can expose environment variables, labels, the definition file, tests, and
runscripts.

```
# Returns a large data structure with all apps
$cseq inspect

# Returns for a particular app
$cseq inspect download
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
```
1. Download data, map the data base to an empty folder on our local machine
     (we actually will do this from the browser as the sri toolkit is changed.
2. Perform mapping step of pipeline, mapping the same folder.
3. perform poisson regression on filtered reads
4. Finally, sort results

Let's add a volume to our command:

```
cseq="docker run -v $PWD/data:/scif/data vanessa/cseq"
```

## CarrierSeq Pipeline
The common user might want access to the larger scoped pipeline that the software provides. In the case of CarrierSeq, this means (optionally, download) mapping, poisson, and then sorting. If the image isn't provided (e.g., a Singularity Registry or Singularity Hub) the user can build from the build recipe file, `Singularity`. If you haven't done this already:

```
docker build -t vanessa/cseq .
```

The author is still working on updating the automated download, for now download from [here](https://www.dropbox.com/sh/vyor82ulzh7n9ke/AAC4W8rMe4z5hdb7j4QhF_IYa?dl=0) and then move into some data folder you intend to mount:

``
mv $HOME/Downloads/all_reads.fastq data/
```

Let's be lazy and put the bind and singularity command into an environment variable so we don't need
to type it many times.

```
cseq="docker run -v $PWD/data:/scif/data vanessa/cseq"
```

And then the various steps of the pipeline are run as was specified above:

```
$cseq run mapping
$cseq run poisson
$cseq run sorting
```

### 1. Mapping
To run mapping, bind the data folder to the image, and specify the app to be mapping:

```
cseq="docker run -v $PWD/data:/scif/data vanessa/cseq"
$cseq run mapping
```

## 2. Poisson

```
cseq="docker run -v $PWD/data:/scif/data vanessa/cseq"
$cseq run poisson
```

## 3. Sorting

```
cseq="docker run -v $PWD/data:/scif/data vanessa/cseq"
$cseq run sorting
```


## How Can I Change It?
Given two containers that share the same input and output schema, I could
swap out of the steps:


```
...
docker run -v $PWD/data:/scif/data vanessa/cseq run sorting
docker run -v $PWD/data:/scif/data vanessa/another run sorting
```

or even provide a single container with multiple options for the same step


```
...
docker run -v $PWD/data:/scif/data vanessa/cseq run sorting1
docker run -v $PWD/data:/scif/data vanessa/cseq run sorting2
```

As a user, you want a container that exposes enough metadata to run different steps of the pipeline, but you don't want to need to know the specifics of each command call or path within the container. In the above, I can direct the container to my mapped input directory
and specify a step in the pipeline, and I dont need to understand how to use `bwa` or `grep` or `seqtk`, or any of the other software
that makes up each.


### Interactive Environment
If you are curious about all the Scientific Filesystem variables available to you:

```
$CSEQ exec mapping env | grep SCIF_
```

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
