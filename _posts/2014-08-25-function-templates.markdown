---
layout: post
title:  "New feature: @function_template"
date:   2014-08-25 21:00:00
categories: gwf update
---

I'm working on a new feature that lets you write templates as Python functions.

The idea is that, instead of writing templates as shell scripts, you can write a template as a Python function that will then automatically get wrapped in a shell script when you create a target and running the target amounts to running the Python function with specified arguments.

There is nothing in this feature that you cannot achieve by writing separate Python scripts and targets that calls these scripts, but that involves at least writing a separate Python script and calling it in a target, and the function template just does this automatically for you.

A motivating example
--------------------

The feature is motivated by a project I am working on right now, where I need to count the different site patterns for a number of different quartets of bears.

I need to count the site patterns for each scaffold in an assembly for each quartet involving two specific bears paired with all other pairs.

I enumerate the quartets to consider like this:

{% highlight python %}
bears =  ["AK017", "AK034", "BB020", "BB034", "BB037", "BB049",
          "BB059", "CON001", "EB027", "WB039"]

def makeQuartets():
    def makeQuartetsGenerator():
        idx1, idx2 = 0, 1
        for idx3 in xrange(2,len(bears)-1):
            for idx4 in xrange(idx3+1,len(bears)):
                yield (idx1, idx2, idx3, idx4)
    return list(makeQuartetsGenerator())

quartets = makeQuartets()
{% endhighlight %}

and want to map over the quartets and do the counting.

Naturally, I want to write a Python function for doing this counting.

I have the genotype counts in a VCF file and information about the scaffolds in a separate file, so the counting is simple enough matter of running through the scaffolds, extracting each scaffold from the VCF file and count the site patterns.

This is where I use the new feature like this:

{% highlight python %}
@function_template(input=['AK_bears.reduced.raw.snps.indels.all.vcf.gz',
                          'AK_bears.reduced.raw.snps.indels.all.vcf.gz.tbi',
                          'scaffold-lengths.txt'])
def countQuartet(quartet, output_file):

    import vcf
    idx1, idx2, idx3, idx4 = quartet

    def count_patterns(reader):
        counts = {'a|bcd': 0, 'b|acd': 0, 'c|abd': 0, 'd|abc': 0, 
                  'ab|cd': 0, 'ac|bd': 0, 'ad|bc': 0}

        for record in reader:
            a = record.samples[idx1].gt_type
            b = record.samples[idx2].gt_type
            c = record.samples[idx3].gt_type
            d = record.samples[idx4].gt_type

            if None in (a, b, c, d): continue

            counts['a|bcd'] += (2-a)*b*c*d + a*(2-b)*(2-c)*(2-d)
            counts['b|acd'] += a*(2-b)*c*d + (2-a)*b*(2-c)*(2-d)
            counts['c|abd'] += a*b*(2-c)*d + (2-a)*(2-b)*c*(2-d)
            counts['d|abc'] += a*b*c*(2-d) + (2-a)*(2-b)*(2-c)*d

            counts['ab|cd'] += (2-a)*(2-b)*c*d + a*b*(2-c)*(2-d)
            counts['ac|bd'] += (2-a)*b*(2-c)*d + a*(2-b)*c*(2-d)
            counts['ad|bc'] += (2-a)*b*c*(2-d) + a*(2-b)*(2-c)*d

        return counts

    reader = vcf.Reader(open('AK_bears.reduced.raw.snps.indels.all.vcf.gz'))
    out = open(output_file, 'w')
    print >> out, 'scaffold length a|bcd b|acd c|abd d|abc ab|cd ac|bd ad|bc'
    with open('scaffold-lengths.txt') as scaffold_file:
        for line in scaffold_file:
            scaffold, length = line.split()
            try:
                counts = count_patterns(reader.fetch(scaffold))

                print >> out, scaffold, length,
                print >> out, counts['a|bcd'], counts['b|acd'], 
                print >> out, counts['c|abd'], counts['d|abc'],
                print >> out, counts['ab|cd'], counts['ac|bd'], counts['ad|bc']

            except ValueError:
                pass
{% endhighlight %}

The function takes the quartet as an argument, holding the indices of the quartet, together with the output file it should write the results to.

It then does the counting.

The magic happens with the `@function_template` adaptor. This works very much like a `template()` function and takes the same arguments as the `template()` function but adapts a Python function making it into a template.

I can create targets from it like this:

{% highlight python %}
for i, q in enumerate(quartets):
    target('countQuartet-{}'.format(i)) << \
        countQuartet(q, 'counts/{}.txt'.format('-'.join(bears[x] for x in q)))
{% endhighlight %}

Calling the function works like calling a template and providing it to a target using the `<<` operator creates the target.

Evaluating the target just runs the function with the parameters provided here.


Limitations
-----------

Since the feature works by serializing the function and the arguments the closure of the function is lost, so any data it needs to work on has to be provided as an argument and any module it needs to work on has to be imported within the function.

This is why I import the `vcf` module inside the function and why the counting function is nested inside `countQuartet()`. Imported modules outside the function won't work since the modules aren't imported in the generated script for the target and functions defined outside the template function won't be visible either.

Until I figure out how to deal with that, you cannot provide default parameters to the function either.

Default parameters might be a hack around closures but the way I wrap the function right now doesn't handle it.


Try it out
----------

To play with the new feature you need to check out the `feature/function_templates` branch.

I have only tested it on the project I am working on right now, but with some more testing I think it could be a cool new feature.


Future work
-----------

I haven't experimented that much with the serialization so I am sure some of the limitations can be removed.

I plan to use this feature to make it simpler to write workflows using the split-map-reduce design pattern.

That is the pattern I am using in the example above. I split the work I need to do into the quartets, then handle each quartet by mapping the function template over the quartets, and later I will combine the result in a "reduce" step, but I need to do some more processing on them that I haven't implemented yet.

The serialization and wrapping of targets is a little slow right now, but I am not caching anything and write the scripts each time the workflow is evaluated so there are plenty of places where it can be optimized.

As always, I welcome comments and suggestions.
