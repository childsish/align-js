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
var GlobalAligner = (function () {
    function GlobalAligner(scoring_matrix, character_map) {
        if (scoring_matrix.length != Math.pow(character_map.alphabet.length, 2)) {
            throw new Error('Scoring matrix must be square of alphabet length');
        }
        this._scoring_matrix = scoring_matrix;
        this._character_map = character_map ? character_map : new CharacterMap(DEFAULT_NUCLEOTIDE_ALPHABET);
    }
    GlobalAligner.prototype.align = function (sequence1, sequence2) {
        var scores = new Int32Array([0, 0, 0]);
        var alignment = new GlobalAlignment(sequence1, sequence2);
        var s1 = this._character_map.translate(sequence1);
        var s2 = this._character_map.translate(sequence2);
        var gap = this._character_map.translate('_')[0];
        for (var i = 0; i < sequence1.length; ++i) {
            for (var j = 0; j < sequence2.length; ++j) {
                scores[GlobalAligner.DIAG] = alignment.get_score_at(i - 1, j - 1) + this._get_score(s1[i], s2[j]);
                scores[GlobalAligner.LEFT] = alignment.get_score_at(i, j - 1) + this._get_score(gap, s2[j]);
                scores[GlobalAligner.UP] = alignment.get_score_at(i - 1, j) + this._get_score(s1[i], gap);
                var idx = scores.indexOf(Math.max.apply(Math, scores));
                alignment.set_score_at(i, j, scores[idx]);
                alignment.set_pointer(i, j, idx);
            }
        }
        return alignment;
    };
    GlobalAligner.prototype._get_score = function (i, j) {
        return this._scoring_matrix[i * this._character_map.alphabet.length + j];
    };
    GlobalAligner.DIAG = 0;
    GlobalAligner.LEFT = 1;
    GlobalAligner.UP = 2;
    return GlobalAligner;
}());
var SemiGlobalAligner = (function () {
    function SemiGlobalAligner(scoring_matrix, character_map) {
        if (scoring_matrix.length != Math.pow(character_map.alphabet.length, 2)) {
            throw new Error('Scoring matrix must be square of alphabet length');
        }
        this._scoring_matrix = scoring_matrix;
        this._character_map = character_map ? character_map : new CharacterMap(DEFAULT_NUCLEOTIDE_ALPHABET);
    }
    SemiGlobalAligner.prototype.align = function (sequence1, sequence2) {
        var scores = new Int32Array([0, 0, 0]);
        var alignment = new SemiGlobalAlignment(sequence1, sequence2);
        var s1 = this._character_map.translate(sequence1);
        var s2 = this._character_map.translate(sequence2);
        var gap = this._character_map.translate('_')[0];
        for (var i = 0; i < sequence1.length; ++i) {
            for (var j = 0; j < sequence2.length; ++j) {
                scores[GlobalAligner.DIAG] = alignment.get_score_at(i - 1, j - 1) + this._get_score(s1[i], s2[j]);
                scores[GlobalAligner.LEFT] = alignment.get_score_at(i, j - 1) + this._get_score(gap, s2[j]);
                scores[GlobalAligner.UP] = alignment.get_score_at(i - 1, j) + this._get_score(s1[i], gap);
                var idx = scores.indexOf(Math.max.apply(Math, scores));
                alignment.set_score_at(i, j, scores[idx]);
                alignment.set_pointer(i, j, idx);
            }
        }
        return alignment;
    };
    SemiGlobalAligner.prototype._get_score = function (i, j) {
        return this._scoring_matrix[i * this._character_map.alphabet.length + j];
    };
    SemiGlobalAligner.DIAG = 0;
    SemiGlobalAligner.LEFT = 1;
    SemiGlobalAligner.UP = 2;
    return SemiGlobalAligner;
}());
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
        for (var i = 0; i < sequence1.length; ++i) {
            for (var j = 0; j < sequence2.length; ++j) {
                scores[LocalAligner.DIAG] = alignment.get_score_at(i - 1, j - 1) + this._get_score(s1[i], s2[j]);
                scores[LocalAligner.LEFT] = alignment.get_score_at(i, j - 1) + this._get_score(gap, s2[j]);
                scores[LocalAligner.UP] = alignment.get_score_at(i - 1, j) + this._get_score(s1[i], gap);
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
    LocalAligner.END = 0;
    LocalAligner.DIAG = 1;
    LocalAligner.LEFT = 2;
    LocalAligner.UP = 3;
    return LocalAligner;
}());
exports.LocalAligner = LocalAligner;
var Alignment = (function () {
    function Alignment(s1, s2) {
        this._scores = new Int32Array((s1.length + 1) * (s2.length + 1));
        this._pointers = new Int32Array((s1.length + 1) * (s2.length + 1));
        this._stop = [s1.length - 1, s2.length - 1];
        this._s1 = s1;
        this._s2 = s2;
    }
    Alignment.prototype.to_string = function () {
        throw new Error('Invalid use of Alignment class. Use GlobalAlignment, LocalAlignment or SemiGlobalAlignment instead.');
    };
    Alignment.prototype.get_score = function () {
        var _a = this._stop, i = _a[0], j = _a[1];
        return this._scores[this._get_index(i, j)];
    };
    Alignment.prototype.get_score_at = function (i, j) {
        return this._scores[this._get_index(i, j)];
    };
    Alignment.prototype.set_score_at = function (i, j, score) {
        this._scores[this._get_index(i, j)] = score;
    };
    Alignment.prototype.get_pointer = function (i, j) {
        return this._pointers[this._get_index(i, j)];
    };
    Alignment.prototype.set_pointer = function (i, j, pointer) {
        this._pointers[this._get_index(i, j)] = pointer;
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
            this.set_score_at(i - 1, -1, -i);
            this.set_pointer(i - 1, -1, GlobalAligner.UP);
        }
        for (var j = 1; j < s2.length + 1; ++j) {
            this.set_score_at(-1, j - 1, -j);
            this.set_pointer(-1, j - 1, GlobalAligner.LEFT);
        }
    }
    GlobalAlignment.prototype.to_string = function () {
        var si = this._s1;
        var sj = this._s2;
        var _a = this._stop, i = _a[0], j = _a[1];
        var a = [];
        var ai = [];
        var aj = [];
        var pointer = this.get_pointer(i, j);
        while (i >= 0 || j >= 0) {
            if (pointer == GlobalAligner.LEFT) {
                ai.push('-');
                aj.push(sj[j]);
                a.push(' ');
                j -= 1;
            }
            else if (pointer == GlobalAligner.UP) {
                ai.push(si[i]);
                aj.push('-');
                a.push(' ');
                i -= 1;
            }
            else if (pointer == GlobalAligner.DIAG) {
                ai.push(si[i]);
                aj.push(sj[j]);
                a.push(si[i] == sj[j] ? '|' : '.');
                i -= 1;
                j -= 1;
            }
            pointer = this.get_pointer(i, j);
        }
        return ai.reverse().join('') + '\n' + a.reverse().join('') + '\n' + aj.reverse().join('');
    };
    return GlobalAlignment;
}(Alignment));
exports.GlobalAlignment = GlobalAlignment;
var SemiGlobalAlignment = (function (_super) {
    __extends(SemiGlobalAlignment, _super);
    function SemiGlobalAlignment(s1, s2) {
        _super.call(this, s1, s2);
        if (s1.length < s2.length) {
            for (var i = 1; i < s1.length + 1; ++i) {
                this.set_score_at(i - 1, -1, -i);
                this.set_pointer(i - 1, -1, GlobalAligner.UP);
            }
            for (var j = 1; j < s2.length + 1; ++j) {
                this.set_pointer(-1, j - 1, GlobalAligner.LEFT);
            }
        }
        else {
            for (var i = 1; i < s1.length + 1; ++i) {
                this.set_pointer(i - 1, -1, GlobalAligner.UP);
            }
            for (var j = 1; j < s2.length + 1; ++j) {
                this.set_score_at(-1, j - 1, -j);
                this.set_pointer(-1, j - 1, GlobalAligner.LEFT);
            }
        }
    }
    SemiGlobalAlignment.prototype.to_string = function () {
        var si = this._s1;
        var sj = this._s2;
        var _a = this._stop, i = _a[0], j = _a[1];
        var a = [];
        var ai = [];
        var aj = [];
        if (si.length < sj.length) {
            aj.push(sj.substring(j + 1));
        }
        else {
            ai.push(si.substring(i + 1));
        }
        var pointer = this.get_pointer(i, j);
        while (i >= 0 && j >= 0) {
            if (pointer == GlobalAligner.LEFT) {
                ai.push('-');
                aj.push(sj[j]);
                a.push(' ');
                j -= 1;
            }
            else if (pointer == GlobalAligner.UP) {
                ai.push(si[i]);
                aj.push('-');
                a.push(' ');
                i -= 1;
            }
            else if (pointer == GlobalAligner.DIAG) {
                ai.push(si[i]);
                aj.push(sj[j]);
                a.push(si[i] == sj[j] ? '|' : '.');
                i -= 1;
                j -= 1;
            }
            pointer = this.get_pointer(i, j);
        }
        if (i > 0) {
            ai.push(si.substring(0, i + 1));
            aj.push(' '.repeat(i + 1));
            a.push(' '.repeat(i + 1));
        }
        else if (j > 0) {
            ai.push(' '.repeat(j + 1));
            aj.push(sj.substring(0, j + 1));
            a.push(' '.repeat(j + 1));
        }
        return ai.reverse().join('') + '\n' + a.reverse().join('') + '\n' + aj.reverse().join('');
    };
    SemiGlobalAlignment.prototype.set_score_at = function (i, j, score) {
        _super.prototype.set_score_at.call(this, i, j, score);
        var on_short_sequence_border = (this._s1.length < this._s2.length && i == this._s1.length - 1) ||
            (this._s2.length < this._s1.length && j == this._s2.length - 1);
        if (on_short_sequence_border && score > this.get_score_at(this._stop[0], this._stop[1])) {
            this._stop = [i, j];
        }
    };
    return SemiGlobalAlignment;
}(Alignment));
var LocalAlignment = (function (_super) {
    __extends(LocalAlignment, _super);
    function LocalAlignment(s1, s2) {
        _super.call(this, s1, s2);
    }
    LocalAlignment.prototype.to_string = function () {
        var si = this._s1;
        var sj = this._s2;
        var _a = this._stop, i = _a[0], j = _a[1];
        var a = [];
        var ai = [];
        var aj = [];
        var pointer = this.get_pointer(i, j);
        while (pointer != LocalAligner.END) {
            if (pointer == LocalAligner.LEFT) {
                ai.push('-');
                aj.push(sj[j]);
                a.push(' ');
                j -= 1;
            }
            else if (pointer == LocalAligner.UP) {
                ai.push(si[i]);
                aj.push('-');
                a.push(' ');
                i -= 1;
            }
            else if (pointer == LocalAligner.DIAG) {
                ai.push(si[i]);
                aj.push(sj[j]);
                a.push(si[i] == sj[j] ? '|' : '.');
                i -= 1;
                j -= 1;
            }
            pointer = this.get_pointer(i, j);
        }
        return ai.reverse().join('') + '\n' + a.reverse().join('') + '\n' + aj.reverse().join('');
    };
    LocalAlignment.prototype.set_score_at = function (i, j, score) {
        _super.prototype.set_score_at.call(this, i, j, score);
        if (score > this.get_score_at(this._stop[0], this._stop[1])) {
            this._stop = [i, j];
        }
    };
    return LocalAlignment;
}(Alignment));
exports.LocalAlignment = LocalAlignment;
function global(sequence1, sequence2) {
    return (new GlobalAligner(DEFAULT_NUCLEOTIDE_SCORING_MATRIX, new CharacterMap(DEFAULT_NUCLEOTIDE_ALPHABET))).align(sequence1, sequence2);
}
exports.global = global;
function semi(sequence1, sequence2) {
    return (new SemiGlobalAligner(DEFAULT_NUCLEOTIDE_SCORING_MATRIX, new CharacterMap(DEFAULT_NUCLEOTIDE_ALPHABET))).align(sequence1, sequence2);
}
exports.semi = semi;
function local(sequence1, sequence2) {
    return (new LocalAligner(DEFAULT_NUCLEOTIDE_SCORING_MATRIX, new CharacterMap(DEFAULT_NUCLEOTIDE_ALPHABET))).align(sequence1, sequence2);
}
exports.local = local;
//# sourceMappingURL=align.js.map