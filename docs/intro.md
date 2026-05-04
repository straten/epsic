# Introduction

epsic is both a library of C++ code that can be used to simulate the 
polarization of electromagnetic radiation and an application that configures
and executes a simulation as described by its command-line arguments.

It was written by Willem van Straten
([\@straten](https://github.com/straten){.user-mention}) in consultation
with Caterina Tiburzi to support our paper, [The Statistics of Radio
Astronomical Polarimetry: Disjoint, Superposed, and Composite
Samples](http://dx.doi.org/10.3847/1538-4357/835/2/293), published in
*The Astrophysical Journal*, 835:293 (2017).

## Disjoint, Superposed, and Composite Samples

epsic can simulate the effects of integration over finite
samples when more than one source is present. The sources can be
disjoint (mutually exclusive), such that only one source emits at a
given instant, or superposed, such that the electric fields are summed.
It is also possible to simulate integration over a composite sample of
unresolved disjoint modes.

## Amplitude modulation

By default, the components of the electric field vector will be normally
distributed; optionally, the amplitude of the vector can be modulated
with an independent random variable drawn from a lognormal distribution.
This amplitude modulating function can be boxcar smoothed or it can be
set to a contiguous sequence of rectangular pulses.

## Overview

epsic performs the following Monte Carlo simulation.

1.  Generate a sequence of $M$ random electric field vector instances
    $\boldsymbol{e}$, each with statistically independent and
    identically distributed (iid) circular complex normal components;
    such a sequence is described by the population mean Stokes
    parameters \[1,0,0,0\].

2.  To yield the desired population mean Stokes parameters, $S_\mu$,
    transform each electric field vector instance by the Hermitian
    square root of
    $2\boldsymbol{\rho}=S_\mu\,\boldsymbol{\sigma}_{\mu}$.

3.  Optionally perform amplitude modulation by multiplying each instance
    of $\boldsymbol{e}$ by an iid random variate $u$ that is drawn
    from a log-normal distribution. The log-normally distributed variate
    is generated from a normally distributed iid variate with zero mean
    and standard deviation $\varsigma$ and is normalized by the mean of
    the distribution, $\langle u \rangle = \exp(\varsigma^2/2)$, 
    such that the mean of the amplitude modulating function is unity.

    -   To simulate rectangular subpulses defined by the sub-sample size
        $n$, a single value of $u$ is applied to $n$ consecutive
        instances of $\boldsymbol{e}$.

    -   The modulation function can also be smoothed using a box-car of
        specified width.

4.  If simulating *superposed samples*, repeat all of the previous steps
    to produce $M$ instances of the electric field vector in the other
    mode then, for each instance of the electric field vectors from
    modes $A$ and $B$, produce $M$ new instances
    $\boldsymbol{e}=\boldsymbol{e}_A+\boldsymbol{e}_B$.

5.  Compute the instantaneous Stokes parameters,
    $s_\mu=\boldsymbol{e}^\dagger \boldsymbol{\sigma}_{\mu} \boldsymbol{e}$.

6.  Optionally divide the sequence of $M$ instantaneous Stokes vectors
    into mutually exclusive Stokes samples of $n$ instances, yielding a
    sequence of $N=M/n$ Stokes samples. This step is not optional when
    simulating composite samples.

7.  If simulating *composite samples*, replace $(1-f)n$ instances in
    each Stokes sample with instantaneous Stokes vectors in the other
    mode.

8.  If simulating *disjoint samples*, replace $(1-F)N$ Stokes samples
    with Stokes samples that contain only instantaneous Stokes vectors
    in the other mode.

9.  For each Stokes sample, compute the sample mean Stokes parameters
    $\bar{S}_\mu$.

10. Compute the $4\times4$ covariances between the Stokes parameters
    using either the $M$ instantaneous Stokes parameters or the $N$
    sample mean Stokes parameters.

Following the simulation, epsic reports

-   the expected population mean Stokes parameters and $4\times4$ matrix
    of covariances between the instantaneous or sample mean Stokes
    parameters;

-   the values of the above quantities as derived from the simulated
    data;

-   the mean value of the degree of polarization, derived from the
    simulated instantaneous or sample mean Stokes parameters; and

-   the modulation index derived from simulated data.

epsic currently does not verify that the computed covariance matrix
matches the theoretical prediction within the uncertainty due to noise.

## Simulation Model Parameters 

The simulation can be configured using the following parameters, each of
which corresponds to a command line option. For a brief list of options
and their arguments, run `epsic -h`.

-   By default, epsic simulates the polarized electromagnetic radiation
    of a single source. To simulate a combination of two sources, use
    one of the following three options:

    -   Enable **superposed modes** using the `-S` option, which has no
        arguments.

    -   Enable **disjoint modes** using the `-D` option; the argument to
        this option is the fraction of Stokes samples that occur in mode
        A in the simulated population.

    -   Enable **composite samples of disjoint modes** using the `-C`
        option; the argument to this option is the fraction of electric
        field instances that occur in mode A in each Stokes sample.

-   The **Stokes sample size** can be set using the `-n` option. The
    argument to this option is the integer number of instances of the
    electric field in each Stokes sample. By default, epsic calculates
    the statistics of the instantaneous Stokes parameters (i.e. $n=1$).

-   The **population size** can be set using the `-N` option. The argument
    to this option is the floating point number of mega Stokes samples
    to simulate, where a mega Stokes sample is $2^{20}$ Stokes samples.
    The total number of instances of electric field vectors that will be
    simulated is the population size times the Stokes sample size.

-   The **population mean Stokes parameters** can be set using the `-s`
    option. The argument to this option is a comma separated list of the
    four population mean Stokes parameters, $I,Q,U,V$; e.g. `epsic -s 5.3,0,0,-1.4`
    will simulate polarized noise with a total intensity
    (Stokes $I$) of 5.3 units and circularly polarized intensity (Stokes
    $V$) of -1.4 units. By default, epsic will simulate unpolarized
    noise with population mean intensity equal to unity.

-   Optional **amplitude modulation** is enabled using the `-l` option.
    The argument to this option is the standard deviation $\varsigma$ of
    the normally distributed variate that is used to generate a
    lognormal distribution. For example, `epsic -l 0.2` will cause epsic
    to multiply each electric field instance by $u=\exp(v/2)$, where $v$
    is normally distributed with mean equal to zero and standard
    deviation equal to 0.2. When amplitude modulation is enabled using
    the `-l` option, the following options apply.

    -   The amplitude modulating function can be set to a contiguous
        sequence of **rectangular pulses** using the `-r` option. The
        argument to this option is the integer number of instances of
        the electric field to be multiplied by a single value of $u$.

    -   The amplitude modulating function can be **boxcar smoothed**
        using the `-b` option. The argument to this option is the width of
        the boxcar expressed as the integer number of instances of the
        electric field that it spans.

If simulating a combination of two sources, either the population mean
Stokes parameters or the modulation properties of the second source can
be specified by preceding the argument to any of `-s`, `-l`, `-r` and/or `-b`
with the letter 'B'; e.g. `epsic -S -l B0.5 -b B4`.

## Cross-covariances between the Stokes parameters 

epsic can also report the measured and predicted cross-covariances
between the Stokes parameters as a function of lag. Currently, this
works only for the instantaneous Stokes parameters. To try it out, add
`-X` to the command line; the results will be printed to a text file named
`acf.txt`.

