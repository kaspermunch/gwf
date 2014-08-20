---
layout: page
title: Documentation
permalink: /documentation/
---

Using GWF involves specifying a workflow in one or more files and using the `gwf` script to submit jobs in the workflow to the computer cluser.

Workflows are specified in Python using a special function, `target()`, that specifies a job to be run. Targets depend on zero or more input files and produce zero or more output files. If any input file is younger than any output file then the target is assumed to be out of date and `gwf` will recognize this.

When you run `gwf` in a directory with a workflow specified in a file named `workflow.py` it will figure out which targets needs to be run and schedule them on the cluster. If there are dependencies between the jobs it will schedule them in such a way that a target is only run after all the targets that it depends on have run. If a target is already submitted to the cluster when you run `gwf` it will not be submitted again, but any other job that depends on it will be set to wait for its completion before it can run.

At any time, you can use `gwf --status` to see the status of your workflow. This will show how many targets are up to date, how many are currently running or queued on the cluster, and how many targets still needs to be submitted.

By default, `gwf --status` will only show the status of terminal targets, i.e. targets that no other target depends on. To see other targets, simply run `gwf --status TARGET_NAME`.

Run `gwf --help` to get a list of options for `gwf`.

For ways of configuring GWF run `gwf-config --help`.


Specifying workflows
====================

Workflows are specified in Python. Having a full blown programming language to specify workflows makes GWF very powerful, but simple workflows can still be specified with very simple syntax.

Below is a very simple workflow with a single target that just unzips a file.

{% highlight python %}
from gwf import *

target('UnZipGenome', input='ponAbe2.fa.gz', output='ponAbe2.fa') << '''
gzcat ponAbe2.fa.gz > ponAbe2.fa
'''
{% endhighlight %}

The first line is necessary to import the `target()` function. 

The call to the `target()` function takes three parameters. The first is always the name of the target -- which should be unique for the workflow -- and any other parameters should be named parameters. These parameters will specify the input and output files for the target and any options that should be passed on to the cluster such as the necessary resources for executing the target.

In the example above, only the input and output files are specified. There is one input file, `ponAbe2.fa.gz`, and one output file `ponAbe2.fa`. When more input or output files are needed these are simply specified as a Python list.

The actual execution for the target is specified outside the call to the function, in the string that follows the `<<` operator. The string that follows the operator should contain the shell commands that the target should execute. In the example there is a single command that uses `gzcat` to unzip a FASTA file.

To string together a number of tasks, a workflow simply needs to specify more than one target. The dependencies between targets are then determined by the input and output files. For a slightly more complex example, see below:

{% highlight python %}
from gwf import *

target('UnZipGenome', input='ponAbe2.fa.gz', output='ponAbe2.fa') << '''
gzcat ponAbe2.fa.gz > ponAbe2.fa
'''

target('IndexGenome', input='ponAbe2.fa',
       output=['ponAbe2.amb', 'ponAbe2.ann', 'ponAbe2.pac']) << '''
bwa index -p ponAbe2 -a bwtsw ponAbe2.fa
'''

target('MapReads', input=['ponAbe2.fa', 
                          'ponAbe2.amb', 'ponAbe2.ann', 'ponAbe2.pac',
                          'Masala_R1.fastq.gz', 'Masala_R2.fastq.gz'],
       output='Masala.sorted.rmdup.bam') << '''

bwa mem -t 16 ponAbe2.fa Masala_R1.fastq.gz Masala_R2.fastq.gz | \
    samtools view -Shb - > /scratch/$GWF_JOBID/unsorted.bam

samtools sort -o /scratch/$GWF_JOBID/unsorted.bam /scratch/$GWF_JOBID/sort | \
    samtools rmdup -s - Masala.sorted.rmdup.bam
'''

{% endhighlight %}

This workflow has three targets: `UnZipGenome`, `IndexGenome`, and `MapReads`. The output file from `UnZipGenome`, `ponAbe2.fa`, is used as the input for `IndexGenome` that in turn produces three index files as output, `ponAbe2.amb`, `ponAbe2.ann`, and `ponAbe2.pac`. Both the FASTA file and the three index files are input files to the `MapReads` target in addition to two FASTQ files. The `MapReads` target executes two shell commands -- any number is allows in a target specification -- to produce the final output.

In this workflow, `UnZipGenome` and `IndexGenome` both have other targets depending on them, so `MapReads` is the only terminal target and `gwf --status` would only show the progress in computing the final output file.


Target options
--------------

Additional named options to the `target()` function can be used to set the resources needed to run the target.

Four options are currently recognized:

* `nodes` -- The number of computer nodes needed for the target. The default is 1.
* `cores` -- The number of cores (per node). The default is 1.
* `memory` -- The necessary RAM. This can be specified either in bytes as an integer or as a string in the usual format like "2g" for 2 gigabytes. The default is "4g".
* `walltime` -- The time needed to finish the job. Specified as a string in the "hh:mm:ss" format. The default is "120:00:00".

In the example above, the `MapReads` target asks `bwa` to use 16 threads for read-mapping so this target should ask for 16 cores:

{% highlight python %}
target('MapReads', input=['ponAbe2.fa', 
                          'ponAbe2.amb', 'ponAbe2.ann', 'ponAbe2.pac',
                          'Masala_R1.fastq.gz', 'Masala_R2.fastq.gz'],
       output='Masala.sorted.rmdup.bam',
       cores=16) << '''

bwa mem -t 16 ponAbe2.fa Masala_R1.fastq.gz Masala_R2.fastq.gz | \
    samtools view -Shb - > /scratch/$GWF_JOBID/unsorted.bam

samtools sort -o /scratch/$GWF_JOBID/unsorted.bam /scratch/$GWF_JOBID/sort | \
    samtools rmdup -s - Masala.sorted.rmdup.bam
'''

{% endhighlight %}

Helper functions in the gwf module
----------------------------------

When importing the `gwf` Python module you also import two helper functions that are useful for working with files: `glob` and `shell`.  The first identifies files using [shell expansion][glob] and returns a Python list of the matched files while the second executes a shell command and returns the result as a list of tokens.

The two commands below will thus both return a list of the files in current directory:

{% highlight python %}
from gwf import *

files = glob("*")
the_same_files = shell("ls")
{% endhighlight %}



Templates
=========

Building every target as in the examples above only works when there are a few targets and when the targets are known at the time the workflow is written. Very often there is a large number of targets and often the exact number of targets depend on the data to be analysed.

When this is the case, we can use more Python programming to create the targets.

If, for example, we have several pairs of FASTQ files and want to map them to a set of BAM files, we can iterate over them and build targets for each BAM file.

{% highlight python %}
from gwf import *

R1files = ['Masala_1_R1.fastq.gz', 'Masala_2_R1.fastq.gz']
R2files = ['Masala_1_R2.fastq.gz', 'Masala_2_R2.fastq.gz']

for i, (r1, r2) in enumerate(zip(R1files)):
  target('MapRead_{}'.format(i+1),
         input=['ponAbe2.fa', 'ponAbe2.amb', 'ponAbe2.ann', 'ponAbe2.pac',
                r1, r2],
         output='Masala_{}.sorted.rmdup.bam'.format(i+1),
         cores=16) << '''

  bwa mem -t 16 ponAbe2.fa {r1} {r2} | \
      samtools view -Shb - > /scratch/$GWF_JOBID/unsorted.bam

  samtools sort -o /scratch/$GWF_JOBID/unsorted.bam /scratch/$GWF_JOBID/sort | \
      samtools rmdup -s - Masala_{number}.sorted.rmdup.bam
  '''.format(r1=r1, r2=r2, number=i+1)
{% endhighlight %}

This looks a lot more complicated -- and using templates as described below it can be much simplified -- but we can break it down a bit to see what is going on.

First, we see the paired-end FASTQ files in two different lists.

When there are many such files, we can build the lists using functions like `glob` or collect them some other way, but here they are just shown explicitly for the example.

Using the Python function `zip` we match them up as pairs and combined with `enumerate` we get each pair together with a number. Here it will be 0 and 1 since those are the numbers `enumerate` will give. So the line

{% highlight python %}
for i, (r1, r2) in enumerate(zip(R1files)):
{% endhighlight %}

will loop over the pairs and give each a number.

The number is assigned to the variable `i` and the two FASTQ files are assigned to the variables `r1` and `r2`.

The number is used to make a unique target name for each pair using the [format][format] string method.

{% highlight python %}
  target('MapRead_{}'.format(i+1),
{% endhighlight %}

We are adding one to make the targets numbered from 1 rather than zero. We didn't have to, but it matches the way we have numbered the input files better.

The two FASTQ files in the pair are used in the list of input files

{% highlight python %}
input=['ponAbe2.fa', 'ponAbe2.amb', 'ponAbe2.ann', 'ponAbe2.pac', r1, r2],
{% endhighlight %}

and the number again to specify the output file

{% highlight python %}
output='Masala_{}.sorted.rmdup.bam'.format(i+1),
{% endhighlight %}

and all three variables are used when specifying the shell commands for the target, again using [format][format].

{% highlight python %}
.format(r1=r1, r2=r2, number=i+1)
{% endhighlight %}

In this case with named parameters to make the subsitutions more explicit.





Templates
---------

**FIXME**

Using functions to make templates
---------------------------------

**FIXME**

Example read-mapping workflow
-----------------------------

{% highlight python %}
from gwf import *
from gwf.bwa import bwa_index, bwa_map
from gwf.samtools import samtools_sort

target('UnZipGenome', input='ponAbe2.fa.gz', output='ponAbe2.fa') << '''
gzcat ponAbe2.fa.gz > ponAbe2.fa
'''
target('IndexGenome') << bwa_index(refGenome='ponAbe2')
target('MapReads')    << bwa_map(refGenome='ponAbe2', 
                                 R1='Masala_R1.fastq.gz', 
                                 R2='Masala_R2.fastq.gz', 
                                 bamfile='Masala.unsorted.bam')
target('SortBAM')     << samtools_sort(name='Masala')
{% endhighlight %}



[glob]: http://en.wikipedia.org/wiki/Glob_(programming)
[format]: https://docs.python.org/2/library/string.html#format-string-syntax
