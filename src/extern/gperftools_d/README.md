# gperftools_d
D bindings for gperftools(Google Performance Tools)

<a href="https://code.dlang.org/packages/gperftools_d" title="Go to gperftools_d"><img src="https://img.shields.io/dub/v/gperftools_d.svg" alt="Dub version"></a>
<a href="https://code.dlang.org/packages/gperftools_d" title="Go to gperftools_d"><img src="https://img.shields.io/dub/dt/gperftools_d.svg" alt="Dub downloads"></a>

# Installation

Install google-perftools and graphviz

```sh
sudo apt-get install google-perftools libgoogle-perftools-dev graphviz
```

Install go and then install google-pprof.

```sh
go get github.com/google/pprof
```

# Performance Profiling

Add `gperftools_d` in `dub.json` as a dependency.

```json
  "dependencies": {
    "gperftools_d": "~>0.1.0"
  }
```

Place the code you want to profile within `ProfilerStart()` and `ProfilerStop()`.

Example: In the `examples/source/app.d` file:

```D
import std.stdio;
import gperftools_d.profiler;

int fib(int x) {
  if(x == 0){
    return 0;
  }
  else if(x == 1){
    return 1;
  }
  else{
    return (fib(x-1) + fib(x-2));
  }
}

void main() {
  ProfilerStart();         // Profiling Starts
  foreach (i; 0 .. 30) {
    writeln(fib(i));
  }
  ProfilerStop();          // Profiling Stops
}

```

To profile:

```sh
dub --compiler=ldc2
CPUPROFILE=/tmp/prof.out <path/to/binary> [binary args]
pprof <path/to/binary> /tmp/prof.out      # -pg-like text output
pprof --gv <path/to/binary> /tmp/prof.out # really cool graphical output
pprof --pdf <path/to/binary> /tmp/prof.out > profile.pdf # dump graphical output in profile.pdf
```

### Example output:

[Profile.pdf](https://github.com/prasunanand/gperftools_d/tree/master/examples/profile.pdf)

# LICENSE

This software is distributed under the [BSD 3-Clause License](LICENSE).

Copyright © 2017, Prasun Anand
