"use strict";
var align_1 = require('../align');
var chai_1 = require('chai');
describe('TestBisect', function () {
    it('sequence 1 contained within sequence 2', function () {
        var sequence1 = 'ggccg';
        var sequence2 = 'ataggccggta';
        var alignment = align_1.align(sequence1, sequence2);
        chai_1.expect(alignment.to_string()).to.equal('   ggccg\n   |||||\nataggccggta');
        chai_1.expect(alignment.get_score()).to.equal(5);
    });
    it('sequence 2 contained within sequence 1', function () {
        var sequence1 = 'ataggccggta';
        var sequence2 = 'ggccg';
        var alignment = align_1.align(sequence1, sequence2);
        chai_1.expect(alignment.to_string()).to.equal('ataggccggta\n   |||||\n   ggccg');
        chai_1.expect(alignment.get_score()).to.equal(5);
    });
    it('sequence 1 overlaps sequence 2 left', function () {
        var sequence1 = 'ataggccg';
        var sequence2 = 'ggccgata';
        var alignment = align_1.align(sequence1, sequence2);
        chai_1.expect(alignment.to_string()).to.equal('ataggccg\n   |||||\n   ggccgata');
        chai_1.expect(alignment.get_score()).to.equal(5);
    });
    it('sequence 1 overlaps sequence 2 right', function () {
        var sequence1 = 'ggccgata';
        var sequence2 = 'ataggccg';
        var alignment = align_1.align(sequence1, sequence2);
        chai_1.expect(alignment.to_string()).to.equal('   ggccgata\n   |||||\nataggccg');
        chai_1.expect(alignment.get_score()).to.equal(5);
    });
});
