function terms = polynomialTerms(dofs, order)
    terms = gen_terms(dofs, order);
end

function termCell = gen_terms(numVars, order, prefix)

    if nargin < 3  % prefix is optional input
        prefix = '';
    end
    termCell = cell(1,1);
    if numVars == 1
        for i = 0:order
            termCell{i+1,1} = [prefix, sprintf('x%d^%d', numVars, i)];
        end
    else
        termIndex = 1;
        for i = 0:order
            termCellNext = gen_terms(numVars-1, order-i, [prefix, sprintf('x%d^%d * ', numVars, i)]);
            numNewTerms = size(termCellNext,1);
            termCell(termIndex:termIndex+numNewTerms-1,1) = termCellNext;
            termIndex = termIndex + numNewTerms;
        end
    end
end
