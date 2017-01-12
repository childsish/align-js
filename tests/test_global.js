"use strict";
var align_1 = require('../align');
var chai_1 = require('chai');
describe('TestGlobalAlignment', function () {
    it('has an exact match', function () {
        var sequence1 = 'atacata';
        var sequence2 = 'atagcgcata';
        var alignment = align_1.global(sequence1, sequence2);
        chai_1.expect(alignment.to_string()).to.equal('ata---cata\n|||   ||||\natagcgcata');
        chai_1.expect(alignment.get_score()).to.equal(20);
    });
});
//# sourceMappingURL=test_global.js.map