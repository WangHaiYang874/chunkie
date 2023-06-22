function lvl = initiate_level(chnkr)
    [p,pw,glw] = prolongation();

    if isa(chnkr, 'chunker')
        h = chnkr.h;
    elseif isa(chnkr,'chunkgraph')
        h = [];
        for chunk = chnkr.echnks
            h = [h;chunk.h];
        end
    end

    lvl = chnk.rcip.level(chnkr.r,chnkr.d,h,chnkr.adj,p,pw,glw);
    lvl = lvl.buildstar();
end

function [p,pw,glw] = prolongation()
    [x,w] = lege.exps(16);
    xfin = [x-1;x+1] / 2;
    wfin = [w;w] / 2;
    p = lege.matrin(16,xfin);
    pw = wfin .* p ./ (w');
    glw = w;
end