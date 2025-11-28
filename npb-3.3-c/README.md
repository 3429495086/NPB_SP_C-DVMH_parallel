# Parallelization of NAS Parallel Benchmarks 3.3 (NPB)

This repository represents parallelization process of NAS Parallel Benchmarks 3.3 as a sequence of transform passes. It also contains `Makefile` which simplify analysis and transformation of NPB. TSAR analyzer is used to perform static analysis and source-to-source transformation. DYNA analyzer allows to perform dynamic analysis.

Let us consider the following example:
```bash
make lu CLASS=S ACTION=inline
```

In this example we perform the `inline` action to perform source-level function inlining for LU benchmark. Note, that we should place appropriate directives into a source code to mark calls which should be inlined.

To specify a size of an input of a benchmark `CLASS` is used. Note, that it should be set even in case of source code analysis to avoid compilation errors.

The following default actions is available:

- `exec` - compile a specified benchmark,
- `clean` - remove intermediate files,
- `anls` - perform static analysis with TSAR,
- `dyna` - merge AST-level representation for files of a specified benchmark (files which prevent merge action will not be analyzed), perform IR-level instrumentation and link a specified benchmark with DYNA, then it is possible to run executable file and to obtain results of dynamic analysis,
- `dyna-sroa` - run SROA (Scalar Replacements Of Aggregates) pass before IR-level instrumentation, this disables dynamic analysis of memory references (variables) which can be promoted to be register references,
- `dyna` - perform IR-level instrumentation and link a specified benchmark with DYNA, then it is possible to run executable file and to obtain results of dynamic analysis,
- `llvm-link` - build LLVM IR for all files of a specified benchmark and link all these *.ll files using `llvm-link` to obtain a single *.ll file for the whole benchmark,
- `anls-link` - analyze *.ll file which is a result of `llvm-link` target,
- `dyna-link` - perform IR-level instrumentation for a single *.ll which is a result of `llvm-link` target,
- `dyna-link-sroa` - run SROA (Scalar Replacements Of Aggregates) pass before IR-level instrumentation, this disables dynamic analysis of memory references (variables) which can be promoted to be register references,
- `wc` - compute the number of lines in a specified benchmark,
- `check` - use TSAR to check sources, for example, the `#pragma spf check nomacro` directive allows to check absence of macros in a specified source range,
- `callgraph` - build callgraph for a specified test, then it can be visualized using GraphViz,
- `inline` - perform source-level inlining using TSAR.

All mentioned actions are specified in `sys/make.tsar` file.

The following variables should be set in file `config/make.def`:

- `${LLC}` - compiler for LLVM IR,
- `${PCHC}` - compiler for Clang pre-compiled headers,
- `${TSAR}` - Traits Static Analyzer (TSAR),
- `${OPT}` - LLVM opt tool,
- `${CC}` - compiler for C sources,
- `${CXX}` - compiler for C++ sources,
- `${WC}` - tool which counts number of lines in a file,
- `${DYNA_LIB}` - path to DYNA analyzer,
- `${CFLAGS}`, `${LLFLAGS}`, `${PCHFLAGS}`, `${CXXFLAGS}` - compiler flags,
- `${CLINKFLAGS}` - link-time flags,
- `${TSAR_FLAGS}` - flags which specify how TSAR should emit pre-compiled headers,
- `${TSAR_LINKFLAGS}` - flags which specify how TSAR should merge pre-compiled headers before analysis,
- some other flags are also available.
