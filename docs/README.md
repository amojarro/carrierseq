# CarrierSeq Scientific Filesystem (SCIF)
The Scientific Filesystem (SCIF) is a host (or container) agnostic specification for a filesystem organization, set of environment variables, and functions to control them to make scientific applications discoverable, and easy to use. For both of these examples, we will use a shared SCIF recipe file, [carrierseq.scif](../carrierseq.scif) and install the SCIF with applications for mapping, poisson, and other helpers into **both** a Docker and Singularity container. The recipe could be installed in other containers / hosts, however here we provide two well-known containers for scientific applications.

 - [CarrierSeq Singularity Scientific Filesystem](singularity.md)
 - [CarrierSeq Docker Scientific Filesystem](docker.scif.md)

For each, we will build a development container (with tools like samtools and bwa exposed as entrypoints) and a pipeline production container (to put those commands together to run a step like "mapping."
