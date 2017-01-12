import { local } from '../align'
import { expect } from 'chai'


describe('TestLocalAlignment', function() {
    it('has an exact match', function() {
        let sequence1 = 'ggccg';
        let sequence2 = 'ataggccggta';

        let alignment = local(sequence1, sequence2);
        expect(alignment.to_string()).to.equal('ggccg\n|||||\nggccg');
        expect(alignment.get_score()).to.equal(25);
    });
});
