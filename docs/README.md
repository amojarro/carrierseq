# CarrierSeq Scientific Filesystem (SCIF)
The Scientific Filesystem (SCIF) is a host (or container) agnostic specification for a filesystem organization, set of environment variables, and functions to control them to make scientific applications discoverable, and easy to use. For both of these examples, we will use a shared SCIF recipe file, [carrierseq.scif](../carrierseq.scif) and install the SCIF with applications for mapping, poisson, and other helpers into **both** a Docker and Singularity container. If we need to make changes, we can just edit the one file and both containers would be updated when built. The recipe could be installed in other containers / hosts, however here we provide two well-known containers for scientific applications.

 - [CarrierSeq Singularity Scientific Filesystem](singularity.md)
 - [CarrierSeq Docker Scientific Filesystem](docker.scif.md)

For each, we will build a development container (with tools like samtools and bwa exposed as entrypoints) and a pipeline production container (to put those commands together to run a step like "mapping." For a quick look at how the same commands are run between Docker and Singularity - below we have a Singularity container called `cseq` and a docker container `vanessa/cseq`

## CarrierSeq Pipeline Container
List apps in container

```
$ ./cseq apps
$ docker run vanessa/cseq apps

  download
      help
   mapping
   poisson
    readme
 reference
   sorting

```

Ask for help for the application called "readme"

```
$ ./cseq help readme
$ docker run vanessa/cseq help readme
```

Print the repository's README.md to the console

```
$ ./cseq run readme
$ docker run vanessa/cseq run readme

#### CarrierSeq
#### About

bioRxiv doi: https://doi.org/10.1101/175281

CarrierSeq is a sequence analysis workflow for low-input nanopore
            sequencing which employs a genomic carrier.

           Github Contributors: Angel Mojarro (@amojarro), 
                                Srinivasa Aditya Bhattaru (@sbhattaru), 
                                Christopher E. Carr (@CarrCE), 
                                and Vanessa Sochat (@vsoch).
 
fastq-filter from: https://github.com/nanoporetech/fastq-filter

[MORE]
```

Inspect the container

```
$ ./cseq inspect
$ docker run vanessa/cseq inspect
```

or a single scientific application

```
$ cseq inspect mapping
$ docker run vanessa/cseq inspect mapping
```

Run three entrypoints in a row, and bind a data directory

```
$ singularity run --bind data:/scif/data cseq run mapping
$ singularity run --bind data:/scif/data cseq run poisson
$ singularity run --bind data:/scif/data cseq run sorting

$ docker run -v $PWD/data:/scif/data vanessa/cseq run mapping
$ docker run -v $PWD/data:/scif/data vanessa/cseq run poisson
$ docker run -v $PWD/data:/scif/data vanessa/cseq run sorting
```

Create an interactive session with an application context

```
$ ./cseq shell mapping
$ docker run -it vanessa/cseq shell mapping
```

## CarrierSeq Development Container

List applications

```
$ docker run vanessa/cseq:dev apps
$ ./cseq-dev apps
       bwa
    fqtrim
      help
    python
     seqtk
sra-toolkit
```

Open interactive python

```
$ ./cseq-dev run python
$ docker run -it vanessa/cseq:dev run python
[python] executing /bin/bash /scif/apps/python/scif/runscript
Python 2.7.9 (default, Jun 29 2016, 13:08:31) 
[GCC 4.9.2] on linux2
Type "help", "copyright", "credits" or "license" for more information.
>>> 
```

Load container with bwa on path

```
$ ./cseq-dev shell bwa
$ docker run -it vanessa/cseq:dev shell bwa
[bwa] executing /bin/bash 
$ which bwa
$ /scif/apps/bwa/bin/bwa
```

Wicked cool!
