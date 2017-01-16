# align-js

Global, semi-global and local alignment implemented in Javascript.

## Usage

The simplest way to align two sequences is to import the `align` function from the `align` module. The three possible types of alignment are `global`, `local` and `semi`.

``` javascript
const align = require('align');

var alignment = align('ataggcgata', 'ggcg', type='global');
console.log(alignment.to_string());
console.log(alignment.get_score());

alignment = align('ataggcgata', 'ggcg', type='local');
console.log(alignment.to_string());
console.log(alignment.get_score());

alignment = align('ataggcgata', 'ggcg', type='semi');
console.log(alignment.to_string());
console.log(alignment.get_score());
```

If you want to perform multiple alignments you can create an Aligner. You can specify the type of alignment when instantiating the Aligner. Again, the three possibilities are 'global', 'local' and 'semi'. The default is 'global'.

``` javascript
const Aligner = require('align').Aligner;

// instantiate a global aligner
var global_aligner = new Aligner();
global_aligner = new Aligner('global');

// instantiate a local aligner
var local_aligner = new Aligner('local');

// instantiate a semi-global aligner
var semi_global_aligner = new Aligner('semi');
```

Custom alphabets and scoring matrices can be created used. You will need to defined your own CharacterMap for a new alphabet.

``` javascript
const align = require('align');

var character_map = new align.CharacterMap('abcde');
var scoring_matrix = new Int32Array([
    5, 4, 3, 2, 1,
    4, 5, 4, 3, 2,
    3, 4, 5, 4, 3,
    2, 3, 4, 5, 4,
    1, 2, 3, 4, 5
]);
var aligner = new align.Aligner('local', scoring_matrix, character_map);
var alignment = aligner.align('abcdef', 'aecbfabefacbfe');
console.log(alignment.to_string());
console.log(alignment.get_score());
```

Examples can be found in the examples directory.
