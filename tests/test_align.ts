import { align } from '../align'
import { expect } from 'chai'


describe('TestBisect', function() {
    it('sequence 1 contained within sequence 2', function() {
        let sequence1 = 'ggccg';
        let sequence2 = 'ataggccggta';

        let alignment = align(sequence1, sequence2);
        expect(alignment.to_string()).to.equal('   ggccg\n   |||||\nataggccggta');
        expect(alignment.get_score()).to.equal(5);
    });
    it('sequence 2 contained within sequence 1', function() {
        let sequence1 = 'ataggccggta';
        let sequence2 = 'ggccg';

        let alignment = align(sequence1, sequence2);
        expect(alignment.to_string()).to.equal('ataggccggta\n   |||||\n   ggccg');
        expect(alignment.get_score()).to.equal(5);
    });
    it('sequence 1 overlaps sequence 2 left', function() {
        let sequence1 = 'ataggccg';
        let sequence2 = 'ggccgata';

        let alignment = align(sequence1, sequence2);
        expect(alignment.to_string()).to.equal('ataggccg\n   |||||\n   ggccgata');
        expect(alignment.get_score()).to.equal(5);
    });
    it('sequence 1 overlaps sequence 2 right', function() {
        let sequence1 = 'ggccgata';
        let sequence2 = 'ataggccg';

        let alignment = align(sequence1, sequence2);
        expect(alignment.to_string()).to.equal('   ggccgata\n   |||||\nataggccg');
        expect(alignment.get_score()).to.equal(5);
    });
});
