function S = sessionCounts()
% Sessions per animal, per region, per task. Source of truth for S2A/S2B/S2C/S3B.
%
% Derived from each decoding script's spec.sessionLoaders (one loader per
% candidate session, so the animal is recoverable) and spec.groupMaps (one row
% per region; entry > 0 is the probe that region used, 0 means not recorded)
% in F1H / F2H / F3H / F4H / F5G _decodeTongue.m. F5G's groups are days, not
% regions; all its sessions are tjM1.
%
% Loader lists are longer than session counts where a candidate session yielded
% no data: 26 listed for 22 sessions (Simple), 33 for 26 (Double), 19 for 18 (VTA).
% An animal not recorded in a region is omitted, never entered as 0.

S.tasks   = {'Simple Reward','Delayed Reward','Double Reward','Learning','VTA'};
S.fields  = {'simpleReward','delayedReward','doubleReward','learning','vta'};
S.regions = {'tjM1','ALM','tjS1'};

S.simpleReward.animals  = {'TD10s','TD9s','TD27','TD26'};
S.simpleReward.tjM1     = [3 4 5 6];
S.simpleReward.ALM      = [4 2 3 3];
S.simpleReward.tjS1     = [];

S.delayedReward.animals = {'TD1','TD4','TD13','TD15','TD8','TD22','TD23'};
S.delayedReward.tjM1    = [3 3 3 4 4 4];     % TD1 TD4 TD13 TD15 TD22 TD23
S.delayedReward.ALM     = [4 2 4 5];         % TD1 TD4 TD22 TD23
S.delayedReward.tjS1    = [5 3 6];           % TD13 TD15 TD8

S.doubleReward.animals  = {'YH2','YH1','TD4f','TD7f','TD24','TD25'};
S.doubleReward.tjM1     = [4 5 4 5];         % YH2 YH1 TD4f TD7f
S.doubleReward.ALM      = [3 2 3 3 4];       % YH2 YH1 TD4f TD24 TD25
S.doubleReward.tjS1     = [];

S.learning.animals      = {'TD3l','TD2l','TD4l','TD5l'};
S.learning.tjM1         = [5 5 5 5];
S.learning.ALM          = [];
S.learning.tjS1         = [];

S.vta.animals           = {'TDv1','TDv4','TDv6','TDv5'};
S.vta.tjM1              = [4 4 6 4];
S.vta.ALM               = [];
S.vta.tjS1              = [];

S.nSessions = [22 35 26 20 18];
S.nAnimals  = [ 4  7  6  4  4];
S.nPairs    = [30 50 33 20 18];   % (session, region) pairs the decoding scripts report

for k = 1:numel(S.fields)
    f = S.fields{k};
    assert(numel(S.(f).animals) == S.nAnimals(k), ...
        '%s lists %d animals, expected %d.', S.tasks{k}, numel(S.(f).animals), S.nAnimals(k));
    p = 0;
    for r = 1:numel(S.regions)
        v = S.(f).(S.regions{r});
        assert(all(v > 0) && numel(v) <= S.nAnimals(k), ...
            '%s %s: zero entry, or more animals than the task has.', S.tasks{k}, S.regions{r});
        p = p + sum(v);
    end
    assert(p == S.nPairs(k), ...
        '%s sums to %d pairs, expected %d.', S.tasks{k}, p, S.nPairs(k));
end
end
