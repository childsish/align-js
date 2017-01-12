"use strict";
var DEFAULT_NUCLEOTIDE_ALPHABET = 'ATGCSWRYKMBVHDN';
exports.DEFAULT_NUCLEOTIDE_ALPHABET = DEFAULT_NUCLEOTIDE_ALPHABET;
var DEFAULT_NUCLEOTIDE_SCORING_MATRIX = new Int32Array([
    5, -4, -4, -4, -4, 1, 1, -4, -4, 1, -4, -1, -1, -1, -2,
    -4, 5, -4, -4, -4, 1, -4, 1, 1, -4, -1, -4, -1, -1, -2,
    -4, -4, 5, -4, 1, -4, 1, -4, 1, -4, -1, -1, -4, -1, -2,
    -4, -4, -4, 5, 1, -4, -4, 1, -4, 1, -1, -1, -1, -4, -2,
    -4, -4, 1, 1, -1, -4, -2, -2, -2, -2, -1, -1, -3, -3, -1,
    1, 1, -4, -4, -4, -1, -2, -2, -2, -2, -3, -3, -1, -1, -1,
    1, -4, 1, -4, -2, -2, -1, -4, -2, -2, -3, -1, -3, -1, -1,
    -4, 1, -4, 1, -2, -2, -4, -1, -2, -2, -1, -3, -1, -3, -1,
    -4, 1, 1, -4, -2, -2, -2, -2, -1, -4, -1, -3, -3, -1, -1,
    1, -4, -4, 1, -2, -2, -2, -2, -4, -1, -3, -1, -1, -3, -1,
    -4, -1, -1, -1, -1, -3, -3, -1, -1, -3, -1, -2, -2, -2, -1,
    -1, -4, -1, -1, -1, -3, -1, -3, -3, -1, -2, -1, -2, -2, -1,
    -1, -1, -4, -1, -3, -1, -3, -1, -3, -1, -2, -2, -1, -2, -1,
    -1, -1, -1, -4, -3, -1, -1, -3, -1, -3, -2, -2, -2, -1, -1,
    -2, -2, -2, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
]);
exports.DEFAULT_NUCLEOTIDE_SCORING_MATRIX = DEFAULT_NUCLEOTIDE_SCORING_MATRIX;
var CharacterMap = (function () {
    function CharacterMap(alphabet, case_sensitive) {
        if (case_sensitive === void 0) { case_sensitive = false; }
        this._map = new Int32Array(256);
        this._map.fill(alphabet.length);
        for (var i = 0; i < alphabet.length; ++i) {
            if (case_sensitive) {
                this._map[alphabet.charCodeAt(i)] = i;
            }
            else {
                this._map[alphabet.toUpperCase().charCodeAt(i)] = i;
                this._map[alphabet.toLowerCase().charCodeAt(i)] = i;
            }
        }
        this.alphabet = alphabet;
    }
    CharacterMap.prototype.translate = function (sequence) {
        var code = new Int32Array(sequence.length);
        for (var i = 0; i < sequence.length; ++i) {
            code[i] = this._map[sequence.charCodeAt(i)];
        }
        return code;
    };
    return CharacterMap;
}());
exports.CharacterMap = CharacterMap;
var LocalAligner = (function () {
    function LocalAligner(scoring_matrix, character_map) {
        if (scoring_matrix.length != Math.pow(character_map.alphabet.length, 2)) {
            throw new Error('Scoring matrix must be square of alphabet length');
        }
        this._scoring_matrix = scoring_matrix;
        this._character_map = character_map ? character_map : new CharacterMap(DEFAULT_NUCLEOTIDE_ALPHABET);
    }
    LocalAligner.prototype.align = function (sequence1, sequence2) {
        var scores = new Int32Array([0, 0, 0, 0]);
        var alignment = new LocalAlignment(sequence1, sequence2);
        var s1 = this._character_map.translate(sequence1);
        var s2 = this._character_map.translate(sequence2);
        var gap = this._character_map.translate('_')[0];
        for (var i = 1; i < sequence1.length + 1; ++i) {
            for (var j = 1; j < sequence2.length + 1; ++j) {
                scores[1] = alignment.get_score_at(i, j - 1) + this._get_score(gap, s2[j - 1]);
                scores[2] = alignment.get_score_at(i - 1, j) + this._get_score(s1[i - 1], gap);
                scores[3] = alignment.get_score_at(i - 1, j - 1) + this._get_score(s1[i - 1], s2[j - 1]);
                var idx = scores.indexOf(Math.max.apply(Math, scores));
                alignment.set_score_at(i, j, scores[idx]);
                alignment.set_pointer(i, j, idx);
            }
        }
        return alignment;
    };
    LocalAligner.prototype._get_score = function (i, j) {
        return this._scoring_matrix[i * this._character_map.alphabet.length + j];
    };
    return LocalAligner;
}());
exports.LocalAligner = LocalAligner;
var LocalAlignment = (function () {
    function LocalAlignment(s1, s2) {
        this._scores = new Int32Array((s1.length + 1) * (s2.length + 1));
        this._pointers = new Int32Array((s1.length + 1) * (s2.length + 1));
        this._stop = [0, 0];
        this._s1 = s1;
        this._s2 = s2;
    }
    LocalAlignment.prototype.to_string = function () {
        var si = this._s1;
        var sj = this._s2;
        var _a = this._stop, i = _a[0], j = _a[1];
        var a = [];
        var ai = [si.substring(i, si.length)];
        var aj = [sj.substring(j, sj.length)];
        while (this._scores[this._get_index(i, j)] > 0) {
            var pointer = this._pointers[this._get_index(i, j)];
            if (pointer == 1) {
                ai.push('-');
                aj.push(sj[j - 1]);
                a.push(' ');
                j -= 1;
            }
            else if (pointer == 2) {
                ai.push(si[i - 1]);
                aj.push('-');
                a.push(' ');
                i -= 1;
            }
            else if (pointer == 3) {
                ai.push(si[i - 1]);
                aj.push(sj[j - 1]);
                a.push(si[i - 1] == sj[j - 1] ? '|' : '.');
                i -= 1;
                j -= 1;
            }
        }
        ai.push(si.substring(0, i));
        aj.push(sj.substring(0, j));
        if (i > j) {
            a.push(' '.repeat(i));
            aj.push(' '.repeat(i - j));
        }
        else if (j > i) {
            a.push(' '.repeat(j));
            ai.push(' '.repeat(j - i));
        }
        else {
            a.push(' '.repeat(i));
        }
        return ai.reverse().join('') + '\n' + a.reverse().join('') + '\n' + aj.reverse().join('');
    };
    LocalAlignment.prototype.get_score = function () {
        var _a = this._stop, i = _a[0], j = _a[1];
        return this._scores[this._get_index(i, j)];
    };
    LocalAlignment.prototype.get_score_at = function (i, j) {
        return this._scores[this._get_index(i, j)];
    };
    LocalAlignment.prototype.set_score_at = function (i, j, score) {
        this._scores[this._get_index(i, j)] = score;
        if (score > this._scores[this._get_index(this._stop[0], this._stop[1])]) {
            this._stop = [i, j];
        }
    };
    LocalAlignment.prototype.get_pointer = function (i, j) {
        return this._pointers[this._get_index(i, j)];
    };
    LocalAlignment.prototype.set_pointer = function (i, j, pointer) {
        this._pointers[this._get_index(i, j)] = pointer;
    };
    LocalAlignment.prototype._get_index = function (i, j) {
        return i * (this._s2.length + 1) + j + 1;
    };
    return LocalAlignment;
}());
exports.LocalAlignment = LocalAlignment;
function align(sequence1, sequence2) {
    return (new LocalAligner(DEFAULT_NUCLEOTIDE_SCORING_MATRIX, new CharacterMap(DEFAULT_NUCLEOTIDE_ALPHABET))).align(sequence1, sequence2);
}
exports.align = align;
//# sourceMappingURL=align.js.map