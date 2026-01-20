# lisp-jaynes-cumming-sim

This is my attempt at learning more about circuit QED with lisp. Probably going to be a hard project that takes a while but definitely beneficial for my career so need to do

jc-sim/
├── src/
│   ├── math/
│   │   ├── complex.lisp        ; helpers if needed
│   │   ├── linalg.lisp         ; vectors, matrices
│   │   └── expm.lisp           ; matrix exponential
│   │
│   ├── quantum/
│   │   ├── basis.lisp          ; |n,g>, |n,e> indexing
│   │   ├── operators.lisp      ; a, a†, σ±, σz
│   │   ├── hamiltonian.lisp    ; JC Hamiltonian builder
│   │   └── evolution.lisp      ; time evolution
│   │
│   ├── observables.lisp        ; ⟨σz⟩, ⟨n⟩, entropy
│   └── io/
│       └── export.lisp         ; write CSV / NPZ
│
├── experiments/
│   ├── vacuum-rabi.lisp
│   ├── detuning.lisp
│   └── collapse-revival.lisp
│
├── data/
├── plots/
└── python/
    └── plot.py