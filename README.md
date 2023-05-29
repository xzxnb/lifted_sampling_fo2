# Exact Sampler for Universally Quantified FO2 (UFO2) and Extensions

Exact sampler for FO2 with cardinality constraints.
This tool is for sampling instances or combinatorical structures from FO2 theory with cardinality constraints.


## Input format

- Markov Logic Network (MLN) format, also see [Pracmln](http://www.pracmln.org/mln_syntax.html)
- Cardinality constraint: `|P| = k`

   
### Example input file

A `friends-smokes` MLN:

```
person = {C1, C2, C3, C4, C5, C6}
friends(person,person)
smokes(person)

3 smokes(x)
1 !friends(x,y) v !smokes(x) v smokes(y) # i.e., friends(x,y) ^ smokes(x) => smokes(y). NOTE: only support CNF for now
```

2 colored tree:
```
vertex = 9
E(vertex, vertex)
red(vertex)
black(vertex)

!E(x,x).
!E(x,y) v E(y,x).
red(x) v black(x).
!red(x) v !black(x).
(!E(x,y) v !red(x) v !red(y)) ^ (!E(x,y) v !black(x) v !black(y)).
Tree[E]
```

2 red-black tree with exact k red vertices:
```
vertex = 9
E(vertex, vertex)
red(vertex)
black(vertex)

!E(x,x).
!E(x,y) v E(y,x).
red(x) v black(x).
!red(x) v !black(x).
(!E(x,y) v !red(x) v !red(y)) ^ (!E(x,y) v !black(x) v !black(y)).
Tree[E]
|a| = 4
```

More examples are in [models](models/)


### Installation
Install requirements:
```
$ pip install -r requirements.txt
```
Add path to your PYTHONPATH:
```
$ export PYTHONPATH=$(pwd)/sampling_ufo2:$PYTHONPATH
```


### How to use
Run the following command:
```
$ python sampling_ufo2/sampler.py -i models/friendsmoker.mln -k 10
```
Find more arguments: 
```
$ python sampling_ufo2/sampler.py -h
```

## References