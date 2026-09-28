# Cortical Disengagement — figure code

MATLAB code that reproduces the figures and statistics in the paper from the exported session files released on Zenodo.

Step-by-step instructions, including troubleshooting, are in `HOWTO_run_code.docx`.

> **Paper:** _add citation / DOI_
> **Data:** _add Zenodo DOI_

## Requirements

- MATLAB R2019b or newer (Windows, macOS or Linux)
- Statistics and Machine Learning Toolbox
- Signal Processing Toolbox

No other code is needed. The functions the figure scripts share are in `shared/`. The few functions they borrow from the lab's pipeline are byte-identical copies in `shared/pipelineCopies/` (see `SOURCES.txt` there).

## Setup

1. **Get the code:** clone or download this repository.
2. **Download the data** from Zenodo and unzip it. After unzipping you should have:

   ```
   Data/
     R1/            Simple Reward Task          <ANM>_<DATE>_obj.mat, <ANM>_<DATE>_kin.mat
     R14/           Delayed Reward Task
     R16/           Double Reward Task
     VTA/           VTA Reward Task
     Learning/      Simple Reward Task, first five days of training
     C4Stim_R1/     photoinhibition at contact 4, Simple Reward   (+ <ANM>_<DATE>_laser.mat)
     C4Stim_R16/    photoinhibition at contact 4, Double Reward   (+ _laser.mat)
     GCStim/        photoinhibition at the go cue                 (+ _laser.mat)
   Disengagement Times/
     <ANM>_<DATE>_P<probe>/   HMM-GLM fits (CSV), one folder per session and probe
   ```

3. **Tell the code where the data are.** Open `setPaths.m` and fill in the two lines under `EDIT HERE`:

   ```matlab
   dataRoot = 'D:\CorticalDisengagement\Data';                 % or '~/data/CorticalDisengagement/Data'
   hmmRoot  = 'D:\CorticalDisengagement\Disengagement Times';
   ```

   Alternatively, leave both empty and move the `Data` and `Disengagement Times` folders into the repository folder. Both are git-ignored.

4. **Run a script.** Open any figure script and press **Run**, or type its name with its folder as the current folder. Each script finds `setPaths.m` itself, so you don't need to add anything to the MATLAB path by hand.

`setPaths.m` also sets:

- `cacheRoot`: the decoding and spike-rate caches, and the per-task summary files. The default is `cache/`, which can grow to several GB and is safe to delete.
- `outputRoot`: the CSV tables from `utilities/`. The default is `output/`.

## Script ↔ figure map

The script name gives the panel: `F2H_…` makes Fig. 2H, and `FS4B_…` makes Fig. S4B.

| Folder | Scripts | Panels | Data read |
|---|---|---|---|
| `tongue kinematics/` | `F1B`, `F2B`, `F3B`, `F4C`, `FS8AB` `_lickDuration` | 1B, 2B, 3B, 4C, S8A–B | R1, R14, R16, VTA, Learning |
| `single units/` | `F1D`, `F2D`, `F3D` `_rasterPSTH` | 1D, 2D, 3D | R1, R14, R16 |
| `spike rate heatmaps/` | `F1E`, `F2E`, `F3E`, `F4E`, `F5C` `_rateHeatmap` | 1E, 2E, 3E, 4E, 5C | R1, R14, R16, VTA, Learning |
| `spike rate averaged/` | `F1F`, `F2F`, `F3F`, `F4F` `_spikeRate`, `F5DE_spikeRate` | 1F, 2F, 3F, 4F, 5D–E | R1, R14, R16, VTA, Learning |
| `decoding/` | `F1H`, `F2H`, `F3H`, `F4H`, `F5G` `_decodeTongue` | 1H, 2H, 3H, 4H, 5G | R1, R14, R16, VTA, Learning |
| `decoding/` | `FS4B`, `FS4D`, `FS4F`, `FS4H`, `FS4J` `_decodeJaw` | S4B, S4D, S4F, S4H, S4J | R1, R14, R16, VTA, Learning |
| `decoding/` | `F3H_FS4F_crossTask` | 3H / S4F cross-task statistics | summaries written by F1H, F3H, FS4B, FS4F |
| `optogenetics/` | `F1JK_opto`, `F1LM_opto`, `F3JKL_opto` | 1J–K, 1L–M, 3J–L | C4Stim_R1, GCStim, C4Stim_R16 |
| `hmm eng mode/` | `F2L_engagementHeatmap`, `F2M_engagementMode`, `FS6B_engagementPerAnimal` | 2L, 2M, S6B | R14 + Disengagement Times |
| `spontaneous bouts/` | `FS10_spontaneous`, `FS10_disengageTimeR14` | S10 | Learning; R14 + Disengagement Times |
| `session summaries/` | `S2A_percentCompleted`, `S2BC_timelines`, `S3B_sessionsPerRegion`, `insertionPlot` | S2A, S2B–C, S3B, probe insertion maps | none (values in the files; counts in `sessionCounts.m`) |
| `utilities/` | `countUnits`, `listExcellentUnits` | unit counts reported in the Methods | all tasks |

## Notes

- **Units:** neural analyses include units labeled `good` or `excellent`, with a mean rate above 0.01 Hz (go-cue aligned, −2 to 4 s).
- **Run time:** the decoding scripts are the slow ones, taking tens of minutes to hours per task. They cache the per-session data in `cacheRoot`, so a second run is much faster.
- **Cross-task comparisons:** `F3H_decodeTongue` and `FS4F_decodeJaw` compare against the Simple Reward summary file that `F1H_decodeTongue` / `FS4B_decodeJaw` write to `cacheRoot`. Run those first, or run `F3H_FS4F_crossTask`, which runs all four in order.
- **Figures:** each script opens its figures and prints its statistics to the Command Window. Nothing is saved automatically.

## Repository layout

```
setPaths.m               <- the only file you need to edit
shared/                  loaders for the exported files (slimMeta, slimToLegacy, loadSlimSession, loadBehavSession, ...)
shared/pipelineCopies/   unmodified copies of pipeline functions (findTrials, findClusters, mySmooth, linspecer, ...)
<analysis folders>/      one script per figure panel (table above)
utilities/               unit counts
```

## License

_add license_
