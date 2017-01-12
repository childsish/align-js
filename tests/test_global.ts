import { global } from '../align'
import { expect } from 'chai'


describe('TestGlobalAlignment', function() {
    it('has an exact match', function() {
        let sequence1 = 'atacata';
        let sequence2 = 'atagcgcata';

        let alignment = global(sequence1, sequence2);
        expect(alignment.to_string()).to.equal('ata---cata\n|||   ||||\natagcgcata');
        expect(alignment.get_score()).to.equal(20);
    });
});
