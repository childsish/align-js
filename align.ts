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

class Aligner {

    static readonly END = 0;
    static readonly DIAG = 1;
    static readonly LEFT = 2;
    static readonly UP = 3;

    _mode : string;
    _scoring_matrix : Int32Array;
    _character_map : CharacterMap;

    constructor(mode = 'global', scoring_matrix : Int32Array, character_map : CharacterMap) {
        if (scoring_matrix.length != character_map.alphabet.length ** 2) {
            throw new Error('Scoring matrix must be square of alphabet length');
        }
        this._mode = mode;
        this._scoring_matrix = scoring_matrix;
        this._character_map = character_map ? character_map : new CharacterMap(DEFAULT_NUCLEOTIDE_ALPHABET);
    }

    align(sequence1 : string, sequence2 : string) : GlobalAlignment {
        let end = this._mode == 'local' ? 0 : 1 << 31;
        let scores = new Int32Array([end, 0, 0, 0]);

        let alignment = this._mode == 'semi' ? new SemiGlobalAlignment(sequence1, sequence2) :
            this._mode == 'local' ? new LocalAlignment(sequence1, sequence2) :
            new GlobalAlignment(sequence1, sequence2);
        let s1 = this._character_map.translate(sequence1);
        let s2 = this._character_map.translate(sequence2);
        let gap = this._character_map.translate('_')[0];
        for (let i = 0; i < sequence1.length; ++i) {
            for (let j = 0; j < sequence2.length; ++j) {
                scores[Aligner.DIAG] = alignment.get_entry(i - 1, j - 1).score + this._get_score(s1[i], s2[j]);
                scores[Aligner.LEFT] = alignment.get_entry(i, j - 1).score + this._get_score(gap, s2[j]);
                scores[Aligner.UP] = alignment.get_entry(i - 1, j).score + this._get_score(s1[i], gap);
                let idx = scores.indexOf(Math.max(...scores));
                alignment.set_entry(i, j, scores[idx], idx);
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
    _max : number[];
    _s1 : string;
    _s2 : string;

    constructor(s1 : string, s2 : string) {
        this._scores = new Int32Array((s1.length + 1) * (s2.length + 1));
        this._pointers = new Int32Array((s1.length + 1) * (s2.length + 1));
        this._stop = [s1.length - 1, s2.length - 1];
        this._max = [s1.length - 1, s2.length - 1];

        this._s1 = s1;
        this._s2 = s2;
    }

    to_string() : string {
        let si = this._s1;
        let sj = this._s2;
        let [i, j] = this._stop;

        let a = [];
        let ai = [];
        let aj = [];

        let pointer = this.get_entry(i, j).pointer;
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
    }

    get_score() : number {
        let [i, j] = this._max;
        return this._scores[this._get_index(i, j)];
    }

    set_entry(i : number, j : number, score : number, pointer : number) {
        let index = this._get_index(i, j);
        this._scores[index] = score;
        this._pointers[index] = pointer
    }

    get_entry(i : number, j : number) {
        let index = this._get_index(i, j);
        return {
            score: this._scores[index],
            pointer: this._pointers[index]
        }
    }

    _get_index(i : number, j : number) {
        return (i + 1) * (this._s2.length + 1) + (j + 1);
    }
}

class GlobalAlignment extends Alignment {

    constructor(s1 : string, s2 : string) {
        super(s1, s2);
        for (let i = 1; i < s1.length + 1; ++i) {
            this.set_entry(i - 1, -1, -i, Aligner.UP);
        }
        for (let j = 1; j < s2.length + 1; ++j) {
            this.set_entry(-1, j - 1, -j, Aligner.LEFT);
        }
    }
}

class SemiGlobalAlignment extends Alignment {

    constructor(s1 : string, s2 : string) {
        super(s1, s2);
        if (s1.length < s2.length) {
            for (let i = 1; i < s1.length + 1; ++i) {
                this.set_entry(i - 1, -1, -i, Aligner.UP);
            }
            for (let j = 1; j < s2.length + 1; ++j) {
                this.set_entry(-1, j - 1, 0, Aligner.LEFT);
            }
        }
        else {
            for (let i = 1; i < s1.length + 1; ++i) {
                this.set_entry(i - 1, -1, 0, Aligner.UP);
            }
            for (let j = 1; j < s2.length + 1; ++j) {
                this.set_entry(-1, j - 1, -j, Aligner.LEFT);
            }
        }
    }

    set_entry(i : number, j : number, score : number, pointer : number) {
        if (this._s1.length < this._s2.length && i == this._s1.length - 1) {
            if (score > this.get_score()) {
                super.set_entry(i, j, score, pointer);
                this._max = [i, j];
            }
            else {
                super.set_entry(i, j, score, Aligner.LEFT);
            }
        }
        else if (this._s2.length < this._s1.length && j == this._s2.length - 1) {
            if (score > this.get_score()) {
                super.set_entry(i, j, score, pointer);
                this._max = [i, j];
            }
            else {
                super.set_entry(i, j, score, Aligner.UP);
            }
        }
        else {
            super.set_entry(i, j, score, pointer);
        }
    }
}

class LocalAlignment extends Alignment {

    set_entry(i : number, j : number, score : number, pointer : number) {
        super.set_entry(i, j, score, pointer);
        if (score > this.get_score()) {
            this._stop = [i, j];
            this._max = [i, j];
        }
    }
}

function align(sequence1 : string, sequence2 : string, type = 'global') {
    return (new Aligner(type, DEFAULT_NUCLEOTIDE_SCORING_MATRIX,
        new CharacterMap(DEFAULT_NUCLEOTIDE_ALPHABET))).align(sequence1, sequence2);
}

export {
    DEFAULT_NUCLEOTIDE_ALPHABET,
    DEFAULT_NUCLEOTIDE_SCORING_MATRIX,
    Aligner,
    GlobalAlignment,
    LocalAlignment,
    CharacterMap,
    align
}
