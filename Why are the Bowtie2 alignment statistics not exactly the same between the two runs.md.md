“为什么两次 bowtie2 比对统计结果不完全一样”

---

## 现象
zcat output/04_trim_linker/DNaseC_K562_NaOH_R1r1_trimLk3_1.fq.gz | head -8
@E250150528L1C001R0010001006;GT;no_adapter;no_adapter
ATGGTCTCCCAGTACTGTTTGTAGGTACTGTTTTTCTTGGGTTGGGGTGGGGGCGTAGATGT
+
GCIICICIIIGIBGICICCCICGIICGICICCCCCECCIIICCIIGIB?IIIII?BGIC@IC
@E250150528L1C001R0010006313;GG;no_adapter;no_adapter
TTCCCATAAAAACTAGACTGAAGCCTTCTCAGAAACTTGATCCTCAT
+
CC?DIGCGGGGG7C?GG*@??A?@BCCDCIGHDGG0?C?.?FICIG@

zcat output_DNaseC_K562_NaOH_R1r1/04_trim_linker/DNaseC_K562_NaOH_R1r1_trimLk3_1.fq.gz | head -8
@E250150528L1C001R0010001006/;GT;no_adapter;no_adapter
ATGGTCTCCCAGTACTGTTTGTAGGTACTGTTTTTCTTGGGTTGGGGTGGGGGCGTAGATGT
+
GCIICICIIIGIBGICICCCICGIICGICICCCCCECCIIICCIIGIB?IIIII?BGIC@IC
@E250150528L1C001R0010006313/;GG;no_adapter;no_adapter
TTCCCATAAAAACTAGACTGAAGCCTTCTCAGAAACTTGATCCTCAT
+
CC?DIGCGGGGG7C?GG*@??A?@BCCDCIGHDGG0?C?.?FICIG@

同一套参数、同一套参考基因组，两次运行：
```bash
bowtie2 -p 80 --seed 42 -x /ssd/index/bowtie2/hg38XY+mm10XY \
  -U output/04_trim_linker/DNaseC_K562_NaOH_R1r1_trimLk3_1.fq.gz -S output.sam
```
和
```bash
bowtie2 -p 80 --seed 42 -x /ssd/index/bowtie2/hg38XY+mm10XY \
  -U output_DNaseC_K562_NaOH_R1r1/04_trim_linker/DNaseC_K562_NaOH_R1r1_trimLk3_1.fq.gz -S output_1.sam
```

输出统计分别是（总 reads 数完全一样）：

- 第一次：
  - 116069 aligned 0 times  
  - 1880248 aligned exactly 1 time  
  - 926525 aligned >1 times  
- 第二次：
  - 116113 aligned 0 times  
  - 1880268 aligned exactly 1 time  
  - 926461 aligned >1 times  

整体比对率都是 96.03%，只是每一类的计数相差几十条。

---

## 结论：这是 **正常的、由 Bowtie2 设计导致的随机性差异**

即便你设置了相同的 `--seed 42`，如果两次跑用的 **FASTQ 文件在 read ID 上不完全相同**（比如你的例子，多或少一个 `/`，或者其它 header 格式差异），Bowtie2 对于“多重比对的 reads”会产生略微不同的选择，从而导致“0 次 / 1 次 / >1 次”这些统计数字出现你看到的那种细微差别。

关键点：

1. **Bowtie2 对多重比对使用“伪随机”来打破平局**  
   - 当某个 read 在参考基因组上有多个“得分完全一样好”的比对位置时，Bowtie2 不会固定报同一个，而是通过内部的伪随机数生成器，在这些位置中随机选一个。
2. **这个随机数的种子不仅取决于 `--seed`，还取决于 read 本身的信息**  
   - 种子是多个因素共同决定的函数，包括：
     - read 名称（ID）
     - read 序列
     - 质量值
     - 以及命令行的 `--seed`
   - 所以：
     - **同一个 FASTQ 文件（ID、序列、质量都完全相同）+ 同一个 `--seed` → 结果可完全重现**；
     - 但如果 **read ID 有所不同**，即使序列与质量不变、`--seed` 也相同，内部实际用于每条 read 的“综合种子”也会不一样 → 在多重匹配的候选位置中作出不同的随机选择 → 最终统计略有差异。

3. **差异只会发生在“存在多个等价比对位置”的 reads 上**  
   - 如果某条 read 在基因组上只有一个最优位置（或其它位置明显更差），无论 ID 怎样写，比对结果都基本不会变。
   - 有变化的正是那 31% 多重比对的 reads 中的一小部分，所以你看到：
     - 几万乃至上百万的多重比对 reads 里，只有几十条在“0 次 / 1 次 / >1 次”三个类别之间互相挪动了一点点。

从宏观上看，两次运行的 **整体比对性能与结果分布是等价的**，这些数量级为“几十条”的差异可以认为是 Bowtie2 的正常随机波动，不影响下游大多数分析。

---

## 直接回答

> 两个 fq 除了 read ID 不一样，其他完全一样，为何用 Bowtie 做 mapping 时结果不一样？

一句话概括：

> 因为 Bowtie2 在处理有多个同样好比对位置的 reads 时，会用依赖 **read 名称 + 序列 + 质量 + `--seed`** 的随机种子来随机选择报告哪个位置；你改了 read ID，相当于改了随机种子的一部分，所以这些多重比对 reads 里有一小部分被分到了不同的类别，导致两次统计结果略有差别。  

---
reference https://bowtie-bio.sourceforge.net/bowtie2/manual.shtm
Randomness in Bowtie 2
Bowtie 2's search for alignments for a given read is "randomized." That is, when Bowtie 2 encounters a set of equally-good choices, it uses a pseudo-random number to choose. For example, if Bowtie 2 discovers a set of 3 equally-good alignments and wants to decide which to report, it picks a pseudo-random integer 0, 1 or 2 and reports the corresponding alignment. Arbitrary choices can crop up at various points during alignment.
The pseudo-random number generator is re-initialized for every read, and the seed used to initialize it is a function of the read name, nucleotide string, quality string, and the value specified with --seed. If you run the same version of Bowtie 2 on two reads with identical names, nucleotide strings, and quality strings, and if --seed is set the same for both runs, Bowtie 2 will produce the same output; i.e., it will align the read to the same place, even if there are multiple equally good alignments. This is intuitive and desirable in most cases. Most users expect Bowtie to produce the same output when run twice on the same input.
However, when the user specifies the --non-deterministic option, Bowtie 2 will use the current time to re-initialize the pseudo-random number generator. When this is specified, Bowtie 2 might report different alignments for identical reads. This is counter-intuitive for some users, but might be more appropriate in situations where the input consists of many identical reads.
