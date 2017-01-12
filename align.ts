const DEFAULT_NUCLEOTIDE_ALPHABET = 'ATGCSWRYKMBVHDN_';

const DEFAULT_NUCLEOTIDE_SCORING_MATRIX = new Int32Array([
     5, -4, -4, -4, -4,  1,  1, -4, -4,  1, -4, -1, -1, -1, -2, -5,
    -4,  5, -4, -4, -4,  1, -4,  1,  1, -4, -1, -4, -1, -1, -2, -5,
    -4, -4,  5, -4,  1, -4,  1, -4,  1, -4, -1, -1, -4, -1, -2, -5,
    -4, -4, -4,  5,  1, -4, -4,  1, -4,  1, -1, -1, -1, -4, -2, -5,
    -4, -4,  1,  1, -1, -4, -2, -2, -2, -2, -1, -1, -3, -3, -1, -5,
     1,  1, -4, -4, -4, -1, -2, -2, -2, -2, -3, -3, -1, -1, -1, -5,
     1, -4,  1, -4, -2, -2, -1, -4, -2, -2, -3, -1, -3, -1, -1, -5,
    -4,  1, -4,  1, -2, -2, -4, -1, -2, -2, -1, -3, -1, -3, -1, -5,
    -4,  1,  1, -4, -2, -2, -2, -2, -1, -4, -1, -3, -3, -1, -1, -5,
     1, -4, -4,  1, -2, -2, -2, -2, -4, -1, -3, -1, -1, -3, -1, -5,
    -4, -1, -1, -1, -1, -3, -3, -1, -1, -3, -1, -2, -2, -2, -1, -5,
    -1, -4, -1, -1, -1, -3, -1, -3, -3, -1, -2, -1, -2, -2, -1, -5,
    -1, -1, -4, -1, -3, -1, -3, -1, -3, -1, -2, -2, -1, -2, -1, -5,
    -1, -1, -1, -4, -3, -1, -1, -3, -1, -3, -2, -2, -2, -1, -1, -5,
    -2, -2, -2, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -5,
    -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5
]);

class CharacterMap {

    alphabet : string;
    _map : Int32Array;

    constructor(alphabet, case_sensitive = false) {
        this._map = new Int32Array(256);
        this._map.fill(alphabet.length);
        for (let i = 0; i < alphabet.length; ++i) {
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

    translate(sequence : string) : Int32Array {
        let code = new Int32Array(sequence.length);
        for (let i = 0; i < sequence.length; ++i) {
            code[i] = this._map[sequence.charCodeAt(i)];
        }
        return code;
    }
}

class GlobalAligner {

    static readonly DIAG = 0;
    static readonly LEFT = 1;
    static readonly UP = 2;

    _scoring_matrix : Int32Array;
    _character_map : CharacterMap;

    constructor(scoring_matrix : Int32Array, character_map : CharacterMap) {
        if (scoring_matrix.length != character_map.alphabet.length ** 2) {
            throw new Error('Scoring matrix must be square of alphabet length');
        }
        this._scoring_matrix = scoring_matrix;
        this._character_map = character_map ? character_map : new CharacterMap(DEFAULT_NUCLEOTIDE_ALPHABET);
    }

    align(sequence1 : string, sequence2 : string) : GlobalAlignment {
        let scores = new Int32Array([0, 0, 0]);
        let alignment = new GlobalAlignment(sequence1, sequence2);
        let s1 = this._character_map.translate(sequence1);
        let s2 = this._character_map.translate(sequence2);
        let gap = this._character_map.translate('_')[0];
        for (let i = 0; i < sequence1.length; ++i) {
            for (let j = 0; j < sequence2.length; ++j) {
                scores[GlobalAligner.DIAG] = alignment.get_score_at(i - 1, j - 1) + this._get_score(s1[i], s2[j]);
                scores[GlobalAligner.LEFT] = alignment.get_score_at(i, j - 1) + this._get_score(gap, s2[j]);
                scores[GlobalAligner.UP] = alignment.get_score_at(i - 1, j) + this._get_score(s1[i], gap);
                let idx = scores.indexOf(Math.max(...scores));
                alignment.set_score_at(i, j, scores[idx]);
                alignment.set_pointer(i, j, idx);
            }
        }
        return alignment;
    }

    _get_score(i : number, j : number) : number {
        return this._scoring_matrix[i * this._character_map.alphabet.length + j];
    }
}

class SemiGlobalAligner {

    static readonly DIAG = 0;
    static readonly LEFT = 1;
    static readonly UP = 2;

    _scoring_matrix : Int32Array;
    _character_map : CharacterMap;

    constructor(scoring_matrix : Int32Array, character_map : CharacterMap) {
        if (scoring_matrix.length != character_map.alphabet.length ** 2) {
            throw new Error('Scoring matrix must be square of alphabet length');
        }
        this._scoring_matrix = scoring_matrix;
        this._character_map = character_map ? character_map : new CharacterMap(DEFAULT_NUCLEOTIDE_ALPHABET);
    }

    align(sequence1 : string, sequence2 : string) : SemiGlobalAlignment {
        let scores = new Int32Array([0, 0, 0]);
        let alignment = new SemiGlobalAlignment(sequence1, sequence2);
        let s1 = this._character_map.translate(sequence1);
        let s2 = this._character_map.translate(sequence2);
        let gap = this._character_map.translate('_')[0];
        for (let i = 0; i < sequence1.length; ++i) {
            for (let j = 0; j < sequence2.length; ++j) {
                scores[GlobalAligner.DIAG] = alignment.get_score_at(i - 1, j - 1) + this._get_score(s1[i], s2[j]);
                scores[GlobalAligner.LEFT] = alignment.get_score_at(i, j - 1) + this._get_score(gap, s2[j]);
                scores[GlobalAligner.UP] = alignment.get_score_at(i - 1, j) + this._get_score(s1[i], gap);
                let idx = scores.indexOf(Math.max(...scores));
                alignment.set_score_at(i, j, scores[idx]);
                alignment.set_pointer(i, j, idx);
            }
        }
        return alignment;
    }

    _get_score(i : number, j : number) : number {
        return this._scoring_matrix[i * this._character_map.alphabet.length + j];
    }
}

class LocalAligner {

    static readonly END = 0;
    static readonly DIAG = 1;
    static readonly LEFT = 2;
    static readonly UP = 3;

    _scoring_matrix : Int32Array;
    _character_map : CharacterMap;

    constructor(scoring_matrix : Int32Array, character_map : CharacterMap) {
        if (scoring_matrix.length != character_map.alphabet.length ** 2) {
            throw new Error('Scoring matrix must be square of alphabet length');
        }
        this._scoring_matrix = scoring_matrix;
        this._character_map = character_map ? character_map : new CharacterMap(DEFAULT_NUCLEOTIDE_ALPHABET);
    }

    align(sequence1 : string, sequence2 : string) : LocalAlignment {
        let scores = new Int32Array([0, 0, 0, 0]);
        let alignment = new LocalAlignment(sequence1, sequence2);
        let s1 = this._character_map.translate(sequence1);
        let s2 = this._character_map.translate(sequence2);
        let gap = this._character_map.translate('_')[0];
        for (let i = 0; i < sequence1.length; ++i) {
            for (let j = 0; j < sequence2.length; ++j) {
                scores[LocalAligner.DIAG] = alignment.get_score_at(i - 1, j - 1) + this._get_score(s1[i], s2[j]);
                scores[LocalAligner.LEFT] = alignment.get_score_at(i, j - 1) + this._get_score(gap, s2[j]);
                scores[LocalAligner.UP] = alignment.get_score_at(i - 1, j) + this._get_score(s1[i], gap);
                let idx = scores.indexOf(Math.max(...scores));
                alignment.set_score_at(i, j, scores[idx]);
                alignment.set_pointer(i, j, idx);
            }
        }
        return alignment;
    }

    _get_score(i : number, j : number) : number {
        return this._scoring_matrix[i * this._character_map.alphabet.length + j];
    }
}

class Alignment {

    _scores : Int32Array;
    _pointers : Int32Array;
    _stop : number[];
    _s1 : string;
    _s2 : string;

    constructor(s1 : string, s2 : string) {
        this._scores = new Int32Array((s1.length + 1) * (s2.length + 1));
        this._pointers = new Int32Array((s1.length + 1) * (s2.length + 1));
        this._stop = [s1.length - 1, s2.length - 1];

        this._s1 = s1;
        this._s2 = s2;
    }

    to_string() : string {
        throw new Error('Invalid use of Alignment class. Use GlobalAlignment, LocalAlignment or SemiGlobalAlignment instead.');
    }

    get_score() : number {
        let [i, j] = this._stop;
        return this._scores[this._get_index(i, j)];
    }

    get_score_at(i : number, j : number) : number {
        return this._scores[this._get_index(i, j)];
    }

    set_score_at(i : number, j : number, score : number) {
        this._scores[this._get_index(i, j)] = score;
    }

    get_pointer(i : number, j : number) : number {
        return this._pointers[this._get_index(i, j)];
    }

    set_pointer(i : number, j : number, pointer : number) {
        this._pointers[this._get_index(i, j)] = pointer;
    }

    _get_index(i : number, j : number) {
        return (i + 1) * (this._s2.length + 1) + (j + 1);
    }
}

class GlobalAlignment extends Alignment {

    constructor(s1 : string, s2 : string) {
        super(s1, s2);
        for (let i = 1; i < s1.length + 1; ++i) {
            this.set_score_at(i - 1, -1, -i);
            this.set_pointer(i - 1, -1, GlobalAligner.UP);
        }
        for (let j = 1; j < s2.length + 1; ++j) {
            this.set_score_at(-1, j - 1, -j);
            this.set_pointer(-1, j - 1, GlobalAligner.LEFT);
        }
    }

    to_string() : string {
        let si = this._s1;
        let sj = this._s2;
        let [i, j] = this._stop;

        let a = [];
        let ai = [];
        let aj = [];

        let pointer = this.get_pointer(i, j);
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
    }
}

class SemiGlobalAlignment extends Alignment {

    constructor(s1 : string, s2 : string) {
        super(s1, s2);
        if (s1.length < s2.length) {
            for (let i = 1; i < s1.length + 1; ++i) {
                this.set_score_at(i - 1, -1, -i);
                this.set_pointer(i - 1, -1, GlobalAligner.UP);
            }
            for (let j = 1; j < s2.length + 1; ++j) {
                this.set_pointer(-1, j - 1, GlobalAligner.LEFT);
            }
        }
        else {
            for (let i = 1; i < s1.length + 1; ++i) {
                this.set_pointer(i - 1, -1, GlobalAligner.UP);
            }
            for (let j = 1; j < s2.length + 1; ++j) {
                this.set_score_at(-1, j - 1, -j);
                this.set_pointer(-1, j - 1, GlobalAligner.LEFT);
            }
        }
    }

    to_string() : string {
        let si = this._s1;
        let sj = this._s2;
        let [i, j] = this._stop;

        let a = [];
        let ai = [];
        let aj = [];
        if (si.length < sj.length) {
            aj.push(sj.substring(j + 1));
        }
        else {
            ai.push(si.substring(i + 1));
        }


        let pointer = this.get_pointer(i, j);
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
    }

    set_score_at(i : number, j : number, score : number) {
        super.set_score_at(i, j, score);
        let on_short_sequence_border = (this._s1.length < this._s2.length && i == this._s1.length - 1) ||
            (this._s2.length < this._s1.length && j == this._s2.length - 1);
        if (on_short_sequence_border && score > this.get_score_at(this._stop[0], this._stop[1])) {
            this._stop = [i, j];
        }
    }
}

class LocalAlignment extends Alignment {

    constructor(s1 : string, s2 : string) {
        super(s1, s2);
    }

    to_string() : string {
        let si = this._s1;
        let sj = this._s2;
        let [i, j] = this._stop;

        let a = [];
        let ai = [];
        let aj = [];

        let pointer = this.get_pointer(i, j);
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
    }

    set_score_at(i : number, j : number, score : number) {
        super.set_score_at(i, j, score);
        if (score > this.get_score_at(this._stop[0], this._stop[1])) {
            this._stop = [i, j];
        }
    }
}

function global(sequence1 : string, sequence2 : string) : GlobalAlignment {
    return (new GlobalAligner(DEFAULT_NUCLEOTIDE_SCORING_MATRIX,
        new CharacterMap(DEFAULT_NUCLEOTIDE_ALPHABET))).align(sequence1, sequence2);
}

function semi(sequence1 : string, sequence2 : string) : SemiGlobalAlignment {
    return (new SemiGlobalAligner(DEFAULT_NUCLEOTIDE_SCORING_MATRIX,
        new CharacterMap(DEFAULT_NUCLEOTIDE_ALPHABET))).align(sequence1, sequence2);
}

function local(sequence1 : string, sequence2 : string) : LocalAlignment {
    return (new LocalAligner(DEFAULT_NUCLEOTIDE_SCORING_MATRIX,
        new CharacterMap(DEFAULT_NUCLEOTIDE_ALPHABET))).align(sequence1, sequence2);

}

export {
    DEFAULT_NUCLEOTIDE_ALPHABET,
    DEFAULT_NUCLEOTIDE_SCORING_MATRIX,
    LocalAligner,
    GlobalAlignment,
    LocalAlignment,
    CharacterMap,
    global,
    semi,
    local
}
