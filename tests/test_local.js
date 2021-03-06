"use strict";
var align_1 = require('../align');
var chai_1 = require('chai');
describe('TestLocalAlignment', function () {
    it('has an exact match', function () {
        var sequence1 = 'ggccg';
        var sequence2 = 'ataggccggta';
        var alignment = align_1.local(sequence1, sequence2);
        chai_1.expect(alignment.to_string()).to.equal('ggccg\n|||||\nggccg');
        chai_1.expect(alignment.get_score()).to.equal(25);
    });
});
//# sourceMappingURL=test_local.js.map