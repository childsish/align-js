# align-js
A naive implementation of Smith-Waterman in Javascript

## Usage

The simplest way to align two sequences is to import the align function.

``` javascript
const align = require('align').align;

var alignment = align('ataggcgata', 'ggcg');
console.log(alignment.to_string());
console.log(alignment.get_score());
```

If you want to perform multiple alignments you can create a LocalAligner.

``` javascript
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
var aligner = new LocalAligner(scoring_matrix, character_map);
var alignment = aligner.align('abcdef', 'aecbfabefacbfe');
console.log(alignment.to_string());
console.log(alignment.get_score());
```
