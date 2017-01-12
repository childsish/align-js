const DEFAULT_NUCLEOTIDE_ALPHABET = 'ATGCSWRYKMBVHDN';

const DEFAULT_NUCLEOTIDE_SCORING_MATRIX = new Int32Array([
     5, -4, -4, -4, -4,  1,  1, -4, -4,  1, -4, -1, -1, -1, -2,
    -4,  5, -4, -4, -4,  1, -4,  1,  1, -4, -1, -4, -1, -1, -2,
    -4, -4,  5, -4,  1, -4,  1, -4,  1, -4, -1, -1, -4, -1, -2,
    -4, -4, -4,  5,  1, -4, -4,  1, -4,  1, -1, -1, -1, -4, -2,
    -4, -4,  1,  1, -1, -4, -2, -2, -2, -2, -1, -1, -3, -3, -1,
     1,  1, -4, -4, -4, -1, -2, -2, -2, -2, -3, -3, -1, -1, -1,
     1, -4,  1, -4, -2, -2, -1, -4, -2, -2, -3, -1, -3, -1, -1,
    -4,  1, -4,  1, -2, -2, -4, -1, -2, -2, -1, -3, -1, -3, -1,
    -4,  1,  1, -4, -2, -2, -2, -2, -1, -4, -1, -3, -3, -1, -1,
     1, -4, -4,  1, -2, -2, -2, -2, -4, -1, -3, -1, -1, -3, -1,
    -4, -1, -1, -1, -1, -3, -3, -1, -1, -3, -1, -2, -2, -2, -1,
    -1, -4, -1, -1, -1, -3, -1, -3, -3, -1, -2, -1, -2, -2, -1,
    -1, -1, -4, -1, -3, -1, -3, -1, -3, -1, -2, -2, -1, -2, -1,
    -1, -1, -1, -4, -3, -1, -1, -3, -1, -3, -2, -2, -2, -1, -1,
    -2, -2, -2, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
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

class LocalAligner {

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
        for (let i = 1; i < sequence1.length + 1; ++i) {
            for (let j = 1; j < sequence2.length + 1; ++j) {
                scores[1] = alignment.get_score_at(i, j - 1) + this._get_score(gap, s2[j - 1]);
                scores[2] = alignment.get_score_at(i - 1, j) + this._get_score(s1[i - 1], gap);
                scores[3] = alignment.get_score_at(i - 1, j - 1) + this._get_score(s1[i - 1], s2[j - 1]);
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

class LocalAlignment {

    _scores : Int32Array;
    _pointers : Int32Array;
    _stop : number[];
    _s1 : string;
    _s2 : string;

    constructor(s1 : string, s2 : string) {
        this._scores = new Int32Array((s1.length + 1) * (s2.length + 1));
        this._pointers = new Int32Array((s1.length + 1) * (s2.length + 1));
        this._stop = [0, 0];

        this._s1 = s1;
        this._s2 = s2;
    }

    to_string() : string {
        let si = this._s1;
        let sj = this._s2;
        let [i, j] = this._stop;

        let a = [];
        let ai = [si.substring(i, si.length)];
        let aj = [sj.substring(j, sj.length)];

        while (this._scores[this._get_index(i, j)] > 0) {
            let pointer = this._pointers[this._get_index(i, j)];
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
        if (score > this._scores[this._get_index(this._stop[0], this._stop[1])]) {
            this._stop = [i, j];
        }
    }

    get_pointer(i : number, j : number) : number {
        return this._pointers[this._get_index(i, j)];
    }

    set_pointer(i : number, j : number, pointer : number) {
        this._pointers[this._get_index(i, j)] = pointer;
    }

    _get_index(i : number, j : number) {
        return i * (this._s2.length + 1) + j + 1;
    }
}

function align(sequence1 : string, sequence2 : string) : LocalAlignment {
    return (new LocalAligner(DEFAULT_NUCLEOTIDE_SCORING_MATRIX,
        new CharacterMap(DEFAULT_NUCLEOTIDE_ALPHABET))).align(sequence1, sequence2);
}

export {
    DEFAULT_NUCLEOTIDE_ALPHABET,
    DEFAULT_NUCLEOTIDE_SCORING_MATRIX,
    LocalAligner,
    LocalAlignment,
    CharacterMap,
    align
}
