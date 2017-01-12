import { semi } from '../align'
import { expect } from 'chai'


describe('TestSemiGlobalAlignment', function() {
    it('has an exact match', function() {
        let sequence1 = 'ggccg';
        let sequence2 = 'ataggccggata';

        let alignment = semi(sequence1, sequence2);
        expect(alignment.to_string()).to.equal('   ggccg\n   |||||\nataggccggata');
        expect(alignment.get_score()).to.equal(25);
    });
});
