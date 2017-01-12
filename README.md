# align-js

Global, semi-global and local alignment implemented in Javascript.

## Usage

The simplest way to align two sequences is to import an alignment function from the `align` module. The three possible alignment functions are `global`, `semi` and `local`.

``` javascript
// const align = require('align').global;
// const align = require('align').semi;
const align = require('align').local;

var alignment = align('ataggcgata', 'ggcg');
console.log(alignment.to_string());
console.log(alignment.get_score());
```

If you want to perform multiple alignments you can create an Aligner. Again, there are three possibilities; GlobalAligner, SemiGlobalAligner and LocalAligner.

``` javascript
// const GlobalAligner = require('align').GlobalAligner;
// const SemiGlobalAligner = require('align').SemiGlobalAligner;
const LocalAligner = require('align').LocalAligner;

var aligner = new LocalAligner();
var alignment = aligner.align('ataggcgata', 'ggcg');
console.log(alignment.to_string());
alignment = aligner.align('ggcg', 'ataggcgata');
console.log(alignment.to_string());
```

Custom alphabets can be created using the CharacterMap.

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
var aligner = new align.LocalAligner(scoring_matrix, character_map);
var alignment = aligner.align('abcdef', 'aecbfabefacbfe');
console.log(alignment.to_string());
console.log(alignment.get_score());
```

Examples can be found in the examples directory.
