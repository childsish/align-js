"use strict";
var __extends = (this && this.__extends) || function (d, b) {
    for (var p in b) if (b.hasOwnProperty(p)) d[p] = b[p];
    function __() { this.constructor = d; }
    d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
};
var DEFAULT_NUCLEOTIDE_ALPHABET = 'ATGCSWRYKMBVHDN_';
exports.DEFAULT_NUCLEOTIDE_ALPHABET = DEFAULT_NUCLEOTIDE_ALPHABET;
var DEFAULT_NUCLEOTIDE_SCORING_MATRIX = new Int32Array([
    5, -4, -4, -4, -4, 1, 1, -4, -4, 1, -4, -1, -1, -1, -2, -5,
    -4, 5, -4, -4, -4, 1, -4, 1, 1, -4, -1, -4, -1, -1, -2, -5,
    -4, -4, 5, -4, 1, -4, 1, -4, 1, -4, -1, -1, -4, -1, -2, -5,
    -4, -4, -4, 5, 1, -4, -4, 1, -4, 1, -1, -1, -1, -4, -2, -5,
    -4, -4, 1, 1, -1, -4, -2, -2, -2, -2, -1, -1, -3, -3, -1, -5,
    1, 1, -4, -4, -4, -1, -2, -2, -2, -2, -3, -3, -1, -1, -1, -5,
    1, -4, 1, -4, -2, -2, -1, -4, -2, -2, -3, -1, -3, -1, -1, -5,
    -4, 1, -4, 1, -2, -2, -4, -1, -2, -2, -1, -3, -1, -3, -1, -5,
    -4, 1, 1, -4, -2, -2, -2, -2, -1, -4, -1, -3, -3, -1, -1, -5,
    1, -4, -4, 1, -2, -2, -2, -2, -4, -1, -3, -1, -1, -3, -1, -5,
    -4, -1, -1, -1, -1, -3, -3, -1, -1, -3, -1, -2, -2, -2, -1, -5,
    -1, -4, -1, -1, -1, -3, -1, -3, -3, -1, -2, -1, -2, -2, -1, -5,
    -1, -1, -4, -1, -3, -1, -3, -1, -3, -1, -2, -2, -1, -2, -1, -5,
    -1, -1, -1, -4, -3, -1, -1, -3, -1, -3, -2, -2, -2, -1, -1, -5,
    -2, -2, -2, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -5,
    -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5
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
var Aligner = (function () {
    function Aligner(mode, scoring_matrix, character_map) {
        if (mode === void 0) { mode = 'global'; }
        if (scoring_matrix.length != Math.pow(character_map.alphabet.length, 2)) {
            throw new Error('Scoring matrix must be square of alphabet length');
        }
        this._mode = mode;
        this._scoring_matrix = scoring_matrix;
        this._character_map = character_map ? character_map : new CharacterMap(DEFAULT_NUCLEOTIDE_ALPHABET);
    }
    Aligner.prototype.align = function (sequence1, sequence2) {
        var end = this._mode == 'local' ? 0 : 1 << 31;
        var scores = new Int32Array([end, 0, 0, 0]);
        var alignment = this._mode == 'semi' ? new SemiGlobalAlignment(sequence1, sequence2) :
            this._mode == 'local' ? new LocalAlignment(sequence1, sequence2) :
                new GlobalAlignment(sequence1, sequence2);
        var s1 = this._character_map.translate(sequence1);
        var s2 = this._character_map.translate(sequence2);
        var gap = this._character_map.translate('_')[0];
        for (var i = 0; i < sequence1.length; ++i) {
            for (var j = 0; j < sequence2.length; ++j) {
                scores[Aligner.DIAG] = alignment.get_entry(i - 1, j - 1).score + this._get_score(s1[i], s2[j]);
                scores[Aligner.LEFT] = alignment.get_entry(i, j - 1).score + this._get_score(gap, s2[j]);
                scores[Aligner.UP] = alignment.get_entry(i - 1, j).score + this._get_score(s1[i], gap);
                var idx = scores.indexOf(Math.max.apply(Math, scores));
                alignment.set_entry(i, j, scores[idx], idx);
            }
        }
        return alignment;
    };
    Aligner.prototype._get_score = function (i, j) {
        return this._scoring_matrix[i * this._character_map.alphabet.length + j];
    };
    Aligner.END = 0;
    Aligner.DIAG = 1;
    Aligner.LEFT = 2;
    Aligner.UP = 3;
    return Aligner;
}());
exports.Aligner = Aligner;
var Alignment = (function () {
    function Alignment(s1, s2) {
        this._scores = new Int32Array((s1.length + 1) * (s2.length + 1));
        this._pointers = new Int32Array((s1.length + 1) * (s2.length + 1));
        this._stop = [s1.length - 1, s2.length - 1];
        this._max = [s1.length - 1, s2.length - 1];
        this._s1 = s1;
        this._s2 = s2;
    }
    Alignment.prototype.to_string = function () {
        var si = this._s1;
        var sj = this._s2;
        var _a = this._stop, i = _a[0], j = _a[1];
        var a = [];
        var ai = [];
        var aj = [];
        var pointer = this.get_entry(i, j).pointer;
        while (pointer != Aligner.END) {
            if (pointer == Aligner.LEFT) {
                ai.push('-');
                aj.push(sj[j]);
                a.push(' ');
                j -= 1;
            }
            else if (pointer == Aligner.UP) {
                ai.push(si[i]);
                aj.push('-');
                a.push(' ');
                i -= 1;
            }
            else if (pointer == Aligner.DIAG) {
                ai.push(si[i]);
                aj.push(sj[j]);
                a.push(si[i] == sj[j] ? '|' : '.');
                i -= 1;
                j -= 1;
            }
            pointer = this.get_entry(i, j).pointer;
        }
        return ai.reverse().join('') + '\n' + a.reverse().join('') + '\n' + aj.reverse().join('');
    };
    Alignment.prototype.get_score = function () {
        var _a = this._max, i = _a[0], j = _a[1];
        return this._scores[this._get_index(i, j)];
    };
    Alignment.prototype.set_entry = function (i, j, score, pointer) {
        var index = this._get_index(i, j);
        this._scores[index] = score;
        this._pointers[index] = pointer;
    };
    Alignment.prototype.get_entry = function (i, j) {
        var index = this._get_index(i, j);
        return {
            score: this._scores[index],
            pointer: this._pointers[index]
        };
    };
    Alignment.prototype._get_index = function (i, j) {
        return (i + 1) * (this._s2.length + 1) + (j + 1);
    };
    return Alignment;
}());
var GlobalAlignment = (function (_super) {
    __extends(GlobalAlignment, _super);
    function GlobalAlignment(s1, s2) {
        _super.call(this, s1, s2);
        for (var i = 1; i < s1.length + 1; ++i) {
            this.set_entry(i - 1, -1, -i, Aligner.UP);
        }
        for (var j = 1; j < s2.length + 1; ++j) {
            this.set_entry(-1, j - 1, -j, Aligner.LEFT);
        }
    }
    return GlobalAlignment;
}(Alignment));
exports.GlobalAlignment = GlobalAlignment;
var SemiGlobalAlignment = (function (_super) {
    __extends(SemiGlobalAlignment, _super);
    function SemiGlobalAlignment(s1, s2) {
        _super.call(this, s1, s2);
        if (s1.length < s2.length) {
            for (var i = 1; i < s1.length + 1; ++i) {
                this.set_entry(i - 1, -1, -i, Aligner.UP);
            }
            for (var j = 1; j < s2.length + 1; ++j) {
                this.set_entry(-1, j - 1, 0, Aligner.LEFT);
            }
        }
        else {
            for (var i = 1; i < s1.length + 1; ++i) {
                this.set_entry(i - 1, -1, 0, Aligner.UP);
            }
            for (var j = 1; j < s2.length + 1; ++j) {
                this.set_entry(-1, j - 1, -j, Aligner.LEFT);
            }
        }
    }
    SemiGlobalAlignment.prototype.set_entry = function (i, j, score, pointer) {
        if (this._s1.length < this._s2.length && i == this._s1.length - 1) {
            if (score > this.get_score()) {
                _super.prototype.set_entry.call(this, i, j, score, pointer);
                this._max = [i, j];
            }
            else {
                _super.prototype.set_entry.call(this, i, j, score, Aligner.LEFT);
            }
        }
        else if (this._s2.length < this._s1.length && j == this._s2.length - 1) {
            if (score > this.get_score()) {
                _super.prototype.set_entry.call(this, i, j, score, pointer);
                this._max = [i, j];
            }
            else {
                _super.prototype.set_entry.call(this, i, j, score, Aligner.UP);
            }
        }
        else {
            _super.prototype.set_entry.call(this, i, j, score, pointer);
        }
    };
    return SemiGlobalAlignment;
}(Alignment));
var LocalAlignment = (function (_super) {
    __extends(LocalAlignment, _super);
    function LocalAlignment() {
        _super.apply(this, arguments);
    }
    LocalAlignment.prototype.set_entry = function (i, j, score, pointer) {
        _super.prototype.set_entry.call(this, i, j, score, pointer);
        if (score > this.get_score()) {
            this._stop = [i, j];
            this._max = [i, j];
        }
    };
    return LocalAlignment;
}(Alignment));
exports.LocalAlignment = LocalAlignment;
function align(sequence1, sequence2, type) {
    if (type === void 0) { type = 'global'; }
    return (new Aligner(type, DEFAULT_NUCLEOTIDE_SCORING_MATRIX, new CharacterMap(DEFAULT_NUCLEOTIDE_ALPHABET))).align(sequence1, sequence2);
}
exports.align = align;
//# sourceMappingURL=align.js.map