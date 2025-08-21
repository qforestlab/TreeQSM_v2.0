%Andrew Burt - a.burt.12@ucl.ac.uk

function runCylinderModel(cloud_name,model_names,workers)
    addpath(genpath('/home/zcooper/TreeQSM/v2'));
    notxt = strsplit(cloud_name,{'.txt'});    
    P = load(cloud_name);
    str = ['/home/zcooper/TreeQSM/v2/', notxt{1}, '.mat'];
    save(str, 'P');
    distcomp.feature('LocalUseMpiexec',false);
    pool = parcluster('local');
    pool.NumWorkers = workers;
    parpool(pool,pool.NumWorkers);
    parfor i=1:length(model_names)
        info = strsplit(model_names{i},'-');
        dmin0 = str2double(info{2});
        rcov0 = str2double(info{3});
        nmin0 = str2double(info{4});
        dmin = str2double(info{5});
        rcov = str2double(info{6});
        nmin = str2double(info{7});
        lcyl = str2double(info{8});
        noground = str2double(info{9});
        try
            qsm_tree(P,dmin0,rcov0,nmin0,dmin,rcov,nmin,lcyl,true,model_names{i},1,1,1,1);
        catch
        end
    end
    exit;
end
