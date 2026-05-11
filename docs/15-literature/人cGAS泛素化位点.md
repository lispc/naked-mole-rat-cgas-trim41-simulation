# **人源 cGAS 蛋白特异性泛素化修饰网络与赖氨酸位点全景研究报告**

*This is for informational purposes only. For medical advice or diagnosis, consult a professional. (本报告仅供参考。如需医疗建议或诊断，请咨询专业人士。)*

## **1\. 天然免疫受体 cGAS 的生物学基础与翻译后修饰需求**

在高等真核生物的天然免疫系统（Innate Immune System）中，精准识别入侵病原体并迅速启动防御反应是维持宿主生存的核心机制。环鸟苷酸-腺苷酸合成酶（Cyclic GMP-AMP Synthase, cGAS，亦称 MB21D1）作为一种进化上高度保守的胞质双链 DNA（dsDNA）模式识别受体（PRR），在此过程中发挥着无可替代的枢纽作用1。当病毒、细菌感染或细胞自身发生严重应激（如线粒体损伤、微核破裂或基因组不稳定性增加）时，异常游离的 dsDNA 会暴露于细胞质或游离于核质中。cGAS 能够以非序列依赖的方式识别并结合这些 dsDNA，其带正电荷的表面及特异性的锌指（Zinc-finger）结构域与 DNA 的糖-磷酸骨架发生紧密互作4。

这种物理结合会诱导 cGAS 发生显著的别构效应，促使其催化口袋开放，进而利用细胞内丰富的 ATP 和 GTP 合成一种独特的非经典环二核苷酸第二信使——2'3'-环鸟苷酸-腺苷酸（2'3'-cGAMP）7。随后，cGAMP 作为配体结合并激活锚定于内质网膜上的干扰素基因刺激蛋白（STING）。激活后的 STING 发生寡聚化并向高尔基体发生囊泡转运，在此过程中招募 TANK 结合激酶 1（TBK1）和 IκB 激酶（IKK），最终导致干扰素调节因子 3（IRF3）和核因子 κB（NF-κB）的磷酸化与核易位，驱动 I 型干扰素（Type I IFN）及广泛的促炎细胞因子网络的表达4。

然而，免疫反应的双刃剑特性决定了 cGAS 的活性必须受到极度严密的监控。cGAS-STING 通路的钝化会导致宿主对 DNA 病毒和致瘤突变极度易感；相反，其过度活化或持续性开启则会引发强烈的自身免疫性疾病与炎症性组织损伤，如 Aicardi-Goutières 综合征（AGS）、系统性红斑狼疮（SLE）以及多种严重神经退行性疾病11。为了在免疫防御与免疫耐受之间维持精准的生理稳态，细胞演化出了极为庞大且精细的翻译后修饰（Post-Translational Modifications, PTMs）网络，在蛋白质水平上实时调控 cGAS 的催化活性、底物亲和力、亚细胞定位及其半衰期5。在所有已知 PTMs 中，泛素化（Ubiquitination）修饰因其高度的生化复杂性与结构多样性，构成了决定 cGAS 信号传导命运的“分子开关”。

## **2\. 泛素化修饰的拓扑复杂性与“泛素密码”理论**

泛素化是一种高度保守的可逆共价修饰机制，涉及将由 76 个氨基酸组成的泛素（Ubiquitin）多肽单体连接到靶蛋白的特定氨基酸残基上10。这一过程由三级酶促级联反应驱动：泛素激活酶（E1）在消耗 ATP 的条件下激活泛素并形成高能硫酯键；随后泛素被转移至泛素结合酶（E2）；最终，泛素连接酶（E3）识别特异性底物，并催化泛素分子的 C 端甘氨酸与底物靶蛋白上的赖氨酸（Lysine, K）残基的 ε-氨基形成异肽键17。人类基因组编码了超过 600 种 E3 泛素连接酶，这赋予了泛素化修饰极高的底物特异性16。

泛素化修饰的深层调控意义在于其拓扑结构的多样性。泛素分子自身含有七个内部赖氨酸残基（K6、K11、K27、K29、K33、K48 和 K63）以及一个 N 端甲硫氨酸（M1）18。靶蛋白可以被修饰为单泛素化（Monoubiquitination），也可以在这些内部位点上进一步延伸，形成特定连接方式的多聚泛素链（Polyubiquitination）17。这种不同连接方式形成的结构多态性被称为“泛素密码”（Ubiquitin Code），不同的“密码”会被特定的泛素结合域（UBDs）或受体读取，从而引导底物走向截然不同的命运。

在经典的生化范式中，K48 连接的多聚泛素链是引导靶蛋白进入 26S 蛋白酶体进行降解的通用信号；而 K63、K27 以及单泛素化修饰通常不介导降解，而是作为大分子复合物组装的信号脚手架，通过改变底物的空间构象或调节蛋白质-蛋白质相互作用（PPI）来激活信号转导激酶级联7。针对 cGAS 的研究表明，其氨基酸序列中散布着多个可被特异性 E3 酶识别的赖氨酸位点，这些位点上附加的不同类型的泛素链，决定了 cGAS 在抗病毒免疫、自身耐受、核内 DNA 损伤修复以及细胞自噬途径中的最终表现。

## **3\. 人鼠种属差异与 cGAS 序列的赖氨酸编号界定**

在探究 cGAS 特异性泛素化位点时，一个长期的文献阅读障碍在于实验模型中物种来源序列的差异。大多数先驱性的基础免疫学研究利用小鼠模型（*Mus musculus*）或小鼠源性细胞系阐明了 cGAS 的修饰网络，因此大量文献中引用的赖氨酸编号（如 K335、K173 等）均为鼠源 cGAS（m-cGAS）的残基位置22。

然而，在针对人类疾病（如癌症免疫疗法或系统性自身免疫病）的转化医学研究中，明确人源 cGAS（Human cGAS, h-cGAS）的确切位点具有不可替代的药理学意义。人类 cGAS 蛋白由 522 个氨基酸组成，虽然其 NTase 催化核心和 Mab21 DNA 结合结构域在进化上高度保守，但其 N 端无序区以及部分环区存在物种特异性的插入或缺失1。这种序列对比上的偏移导致小鼠与人类 cGAS 在同源保守的赖氨酸修饰位点编号上存在系统性差异。例如，鼠源序列中的 mK335 在人源序列中精确映射为 hK347；mK271 在人源中则往往指代相同的结构保守区域，但也需要通过同源比对加以严格区分22。

本报告通过交叉比对人类 cGAS 晶体结构（如 PDB: 4KM5）、UniProt 序列数据库注释（如 Q8N884）以及最新的质谱修饰组学数据，对直接影响人源 cGAS 功能的具体泛素化赖氨酸位点进行了严格的重新映射和验证，确保所讨论的物理位点在人类病理生理学体系下的准确性。

## **4\. 驱动 cGAS 活化与信号放大的泛素化网络**

在病原体入侵或急性组织应激的初期，机体需要迅速突破 cGAS 的自抑制状态。一部分特定的 E3 泛素连接酶受上游预警信号的招募，通过在 cGAS 的关键赖氨酸残基上附加具有激活效应的泛素标签，重塑其别构形态，从而极大增强对 DNA 的敏感性与 cGAMP 的合成效率。

### **4.1 TRIM56 与 TRIM41 介导的单泛素化修饰与二聚化促进**

三结构域基序（Tripartite Motif, TRIM）蛋白家族是天然免疫调节网络中最重要的 E3 泛素连接酶家族之一。在 DNA 识别级联反应中，TRIM56 被证实是 cGAS 活化的核心正向调节因子22。研究发现，当细胞内感知到异常的 dsDNA 刺激时，TRIM56 会与 cGAS 发生直接的物理结合，并特异性地催化 cGAS 发生单泛素化（Monoubiquitination）21。

通过质谱测序和点突变分析确证，这一修饰在人源 cGAS 中精确发生于 **Lys-347 (K347)** 位点（对应小鼠序列中的 mK335）21。从结构生物学的视角来看，K347 位点地处 cGAS 催化结构域内且毗邻决定受体二聚化的关键界面。未结合 DNA 的游离 cGAS 倾向于以单体形式存在，催化口袋处于关闭状态。当 K347 发生单泛素化后，泛素分子的空间占位效应不仅没有阻断底物结合，反而产生了一种变构推力，显著降低了 cGAS 形成 2:2 异源四聚体（cGAS-DNA 复合体）的能量势垒。这种由 TRIM56 介导的单泛素化不仅提高了 cGAS 对配体 DNA 的结合亲和力，还使得即使在较低浓度的 dsDNA 刺激下，cGAS 也能高效组装成有活性的二聚体构象，大量合成 cGAMP，进而触发强烈的 I 型干扰素表达和抗 DNA 病毒免疫应答21。

与 TRIM56 具有类似效应的还有另一家族成员 TRIM41（又称 RINCK）。大量独立研究确证 TRIM41 同属于 cGAS 稳态的重要调节者，其能够直接结合并对 cGAS 进行单泛素化修饰，协同增强 cGAMP 的合成效率21。**但截至目前，尚无实验研究通过质谱或 K→R 突变扫描鉴定 TRIM41 泛素化 cGAS 的具体位点。** 进一步的机制解析显示，TRIM41 同样介导抗病毒和抗转座子的广谱防御。特别是在限制 LINE-1 (L1) 逆转录转座子跳跃的过程中，定位于细胞核内的 cGAS 会增强逆转录中间产物 ORF2p 与 TRIM41 之间的结合亲和力，通过 TRIM41 对 ORF2p 的降解来维护基因组的稳定性，并在抗衰老机制中发挥深远影响21。

### **4.2 RNF185 介导的 K27 连接多聚泛素化**

除了单泛素化，非经典的 K27 连接多聚泛素化同样构成了 cGAS 激活网络的重要一环。环指蛋白 185（Ring Finger Protein 185, RNF185）被鉴定为首个催化 cGAS 发生此类多聚泛素化的 E3 连接酶21。

机制研究揭示，在病毒感染状态下，RNF185 被迅速招募至 cGAS，并在 cGAS 的 **Lys-137 (K137)** 和 **Lys-384 (K384)** 两个关键赖氨酸残基上，特异性催化形成 K27 连接的多聚泛素链12。这类非降解型泛素链在细胞质中充当了信号转导的分子脚手架。生化测定表明，K137 和 K384 位点的 K27 泛素化可直接提升 cGAS 的酶促动力学速率。当实验性地敲低 RNF185，或者将 cGAS 的 K137 和 K384 位点突变为精氨酸（R）以阻断泛素化反应时，细胞在受到 DNA 病毒攻击后生成的 cGAMP 水平和干扰素产量均呈现断崖式下降6。

值得引起临床高度重视的是，这种激活机制在病理状态下的失调与自身免疫性疾病的发生具有强烈的因果关联。在系统性红斑狼疮（SLE）等自身免疫病患者的生物样本中，RNF185 的转录和翻译水平均呈现显著上调。这种异常的酶丰度导致 cGAS 在极低阈值的内源性 DNA 片段刺激下依然保持高水平的 K27 泛素化状态，从而驱使免疫系统陷入持续的炎症级联激活中12。

### **4.3 ANKIB1 介导的信号小体 K11 泛素化与接头招募**

除了调控 cGAS 蛋白自身的分子内活性，泛素化网络还承担了 cGAS 复合体向外传递信号时的接头功能。近期的前沿发现确立了 K11 连接的泛素链在介导 cGAS-STING 通路信号转导中的核心地位。以往 K11 链一直被认为主要由 APC/C 复合体介导并在细胞周期的有丝分裂期发挥靶向降解作用34，但研究者在《Nature Cell Biology》上发表的结论打破了这一认知局限。

E3 泛素连接酶 ANKIB1 能够被活化的免疫识别复合体选择性募集。在 cGAS 识别 DNA 并形成活化的信号小体（Signalosome）后，ANKIB1 将特异性的 K11 泛素链附着在该信号转导核心组件上9。这一 K11 泛素修饰提供了一个高度亲和的对接平台（Docking platform），特异性地招募自噬相关的接头蛋白 Optineurin (OPTN)。OPTN 进而将激酶 TBK1 牵引至信号复合体中，从而打通了从 cGAS/STING 到 TBK1 甚至 IRF3 的最终磷酸化通路9。在小鼠模型中，缺失 ANKIB1 会极大削弱宿主诱导干扰素的能力并损害其抵御 HSV-1 病毒的屏障9。这一发现将 K11 泛素化确立为免疫级联传导的“第三个泛素密码”，为干预异常的免疫应答提供了一个精确的药理学控制点36。

## **5\. 介导 cGAS 活性抑制与空间阻断的泛素化机制**

正如天然免疫系统需要灵敏的激活开关以对抗病原微生物，它同样依赖强有力的“刹车”系统来防止因内源性微量 DNA 泄漏而引发的过度炎症反应。在这些负反馈调控回路中，泛素化带来的空间位阻效应表现出了惊人的直接性和高效性。

膜相关 RING-CH 结构域蛋白 8（MARCH8，又称 RNF178 或 c-MIR）在此充当了 cGAS 的核心负调节因子。MARCH8 作为一种主要定位于细胞内膜系统的 E3 泛素连接酶，在免疫耐受、病毒包膜糖蛋白识别及内吞体-溶酶体分选网络中承担广泛作用38。机制研究表明，MARCH8 通过其高度保守的 RING-CH 结构域，直接与处于活性状态或游离状态的 cGAS 酶促核心发生物理接触41。

在发生结合后，MARCH8 会在 cGAS 的 **Lys-411 (K411)** 位点上催化形成长链的 K63 连接多聚泛素化41。在大多数生化情境中，K63 泛素链主要用于传递激酶活化信号；然而在 cGAS 的结构背景下，K411 位点紧邻 cGAS 结合双链 DNA 的关键电荷表面。MARCH8 在该点引入庞大的 K63 泛素聚合物后，产生了不可逾越的空间位阻（Steric hindrance）。这种物理屏障直接切断了 cGAS 利用其正电荷区域捕获 dsDNA 的路径，导致 cGAS 的 DNA 结合能力断崖式下降，后续的 cGAMP 生成过程也因此陷于停滞43。

这一负向调控在动物模型水平得到了证实：在 *March8* 基因敲除小鼠模型中，由于失去了针对 cGAS K411 位点的泛素化“刹车”，小鼠的巨噬细胞及成纤维细胞对细胞质 DNA 和 HSV-1 病毒的感染表现出极度敏感和高亢的免疫激活特征，大量分泌抗病毒细胞因子41。由此可见，靶向 K411 的空间占位泛素化是细胞用来微调先天免疫阈值、避免因内源性损伤引发自身炎症的关键适应性机制。

## **6\. 决定 cGAS 蛋白质稳态的靶向降解与去泛素化挽救**

无论激活或抑制，如果不从蛋白质总量上进行控制，持续的刺激最终仍会导致不可控的后果。因此，通过泛素-蛋白酶体系统（UPS）和自噬-溶酶体途径调控 cGAS 的绝对表达丰度和蛋白质半衰期，构成了免疫监控的最后一道防线。

### **6.1 自噬-溶酶体途径：p62 的识别与 K414 位点调控**

在抗病毒反应的中后期，为了有效清除已完成使命的过度活化的 cGAS，细胞启动了由 K48 连接多聚泛素化主导的蛋白降解程序。不同于多数进入 26S 蛋白酶体的可溶性蛋白，cGAS 由于其极易与核酸形成液-液相分离（LLPS）缩合物的特性，其清除过程更多地依赖于选择性自噬（Selective Autophagy）12。

生化分析显示，位于 cGAS 序列后部的 **Lys-414 (K414)** 发生强烈的 K48 连接多聚泛素化修饰后，该修饰不仅削弱了 cGAS 的催化性能，更成为引导 cGAS 走向死亡的“分子标签”。带有 K48 泛素链的 cGAS 被经典的自噬受体 p62（SQSTM1）高效识别并锚定，随后被选择性包裹入自噬小体中，通过与溶酶体融合实现彻底降解21。

有趣的是，病毒入侵初期，宿主必须尽可能保护 cGAS 免于过早降解。此时，去泛素化酶（DUBs）构成了挽救 cGAS 命运的分子屏障。由 I 型干扰素诱导表达的 TRIM14 本身并非活性的泛素连接酶，而是发挥桥接作用。TRIM14 能够特异性招募去泛素化酶 USP14 至 cGAS 复合体。USP14 定向切割 cGAS K414 位点上的 K48 泛素链，破坏了 p62 对 cGAS 的识别接口，成功将其从自噬降解边缘“解救”出来，从而极大延长了 cGAS 在抗病毒免疫建立阶段的有效半衰期22。此外，去泛素化酶 USP15 亦参与调控 cGAS。USP15 能够在切割泛素链稳定 cGAS 蛋白的同时，通过改变其物理化学性质，显著增强 cGAS 在细胞质中进行液-液相分离（LLPS）的能力，进一步优化 cGAS 凝集物（Droplets）内底物的局部浓度，驱动免疫信号的高效扩增21。

### **6.2 核内 cGAS 的细胞周期调控：CRL5-SPSB3 介导的 K427/K428 降解**

过去几年天然免疫领域的最大概念颠覆之一，是认识到 cGAS 并非只定位于细胞质。事实上，大量的 cGAS 被发现存在于细胞核中，甚至在细胞周期的大部分时间里以极高的亲和力锚定在染色质上4。cGAS 通过与其结构上的高亲和力结合域，牢牢扣附在组蛋白 H2A-H2B 异二聚体形成的酸性补丁（Acidic patch）上。这种深度的核小体拴系既封锁了 cGAS 与基因组 DNA 的结合表面，防止其引发致命的自身免疫攻击，又赋予了核内 cGAS 意想不到的新功能——通过物理阻遏同源重组（HR）因子的招募，干预 DNA 损伤修复，进而影响肿瘤的发生与细胞衰老进程48。

对于正处于旺盛增殖周期中的细胞而言，维持核内高浓度的 cGAS 具有严重的基因组毒性隐患。因此，伴随着有丝分裂退出（Mitotic exit），细胞必须启动一套精确的清道夫机制。此时，CRL5-SPSB3 泛素连接酶复合体成为针对核内 cGAS 进行特异性降解的行刑者21。SPSB3 作为该复合物的底物受体蛋白，能够跨越染色质的空间障碍，精准识别处于核小体结合状态下的 cGAS C 末端高度保守的 Asn-Asn (NN) 微型降解子基序48。

在识别配对后，SPSB3 将泛素化核心组件拉近，直接在 cGAS 的 **Lys-427 (K427)** 和 **Lys-428 (K428)** 这对相邻的赖氨酸残基上，连续组装 K48 连接的靶向多聚泛素链51。带有这对重度泛素化标记的核内 cGAS 迅速被分离出染色质，交由细胞核内的蛋白酶体系统摧毁。通过在 K427 和 K428 的这一靶向清理程序，细胞成功维持了核内极低水平的本底 cGAS 浓度，保障了下一次 DNA 复制的安全性并规避了自体 DNA 引起的慢性炎症48。

### **6.3 USP29 与 USP27X：稳固防御体系的组成型去泛素化**

除了针对急性感染相诱导的 USP14 和 TRIM14 系统，在维持 cGAS 的日常生理稳态和长期存活方面，USP29 与 USP27X 两类去泛素化酶同样扮演着举足轻重的角色。

在静息状态或病毒感染的潜伏期，USP29 能够与 cGAS 发生稳定的组成性（Constitutive）物理结合。USP29 表现出极高的 K48 泛素链水解活性，通过持续削减 cGAS 分子表面（特别是在 **hK271** 等具有泛素化敏感性的赖氨酸区域）自发累积的 K48 降解标签，拮抗了泛素-蛋白酶体系统的背景损耗22。体内实验充分证实了该调控轴的重量级意义：在基因敲除小鼠（*Usp29* \-/-）中，cGAS 蛋白遭遇了失控性的降解坍塌，小鼠对 DNA 病毒感染的敏感度显著增加，几乎无法产生足够水平的 I 型干扰素。而在模拟自身免疫性疾病的 *Trex1* \-/- 小鼠模型中，额外敲除 *Usp29* 竟然可以通过加速 cGAS 的彻底清除，逆转并治愈小鼠原本严重的自身免疫与系统性炎症表型46。

USP27X 同样作为 cGAS 复合体的重要守护者被鉴定出。它不仅可以通过直接水解 cGAS 上的 K48 连接多聚泛素链来促进 cGAS 蛋白稳定，还具备调控抗病毒通路另一重要感受器 RIG-I 的泛素化网络的能力21。这反映出细胞内泛素化清除系统具备高度的网络连通性和冗余度，共同维持着天然免疫模式识别受体的寿命下限。

### **6.4 拮抗 RNF185 的去泛素化：JOSD2 调控 K27 链**

去泛素化不仅用于免除降解，还用于阻断过度激活。针对前文所述 RNF185 催化 cGAS（K137 和 K384）发生 K27 连接并极大提升催化活性的机制，机体利用去泛素化酶 JOSD2 实现了特异性的逆向消除。JOSD2 在不改变 cGAS 蛋白半衰期或空间定位的前提下，极其精准地剪切掉 cGAS 表面的 K27 多聚泛素链，使其重新回到低活性的单体/闭合构象56。这种调控在肿瘤微环境中尤为显著：结直肠癌等恶性肿瘤细胞通过上调巨噬细胞中 JOSD2 的表达，使得肿瘤相关巨噬细胞内的 cGAS-STING 通路被强烈压制，直接促使该类巨噬细胞极化为具有免疫抑制特性的 M2 型表型，从而为癌细胞构建了坚固的免疫逃逸屏障56。

## **7\. 泛素化与其他翻译后修饰的空间竞争与串扰网络**

由于泛素化反应严格依赖赖氨酸 ε-氨基的亲核攻击，而赖氨酸残基同时也是小泛素样修饰（SUMOylation）、乙酰化（Acetylation）及甲基化（Methylation）等众多 PTMs 的发生位点，这就不可避免地在 cGAS 的三维表面上引发了极为激烈的“位点争夺战”（Cross-talk）。

### **7.1 SUMO 化修饰对降解型泛素化的占位保护**

小泛素样修饰蛋白（SUMO）与泛素系统存在广泛的相互拮抗。在 DNA 病毒侵入的极早期阶段，免疫系统绝不允许 cGAS 因非特异性的 K48 泛素化而过早降解。为此，E3 SUMO 连接酶 TRIM38 迅速对小鼠 cGAS 的 mK217 和 mK464 两个位点（对应于人类序列的高保守区域 hK231 和 hK479）实施了大量 SUMO 化修饰12。这一快速响应所形成的 SUMO 大分子包裹不仅掩盖了潜在的泛素化识别靶点，而且在空间上屏蔽了带有破坏性的 E3 泛素连接酶的接近，为 cGAS 建立免疫初期的预警赢得了宝贵的时间窗。然而，当急性感染得到控制后，随着时间的推移，特异性蛋白酶 SENP2 会逐渐水解这些保护性的 SUMO 标签，使原本被遮蔽的赖氨酸再次裸露，随后 cGAS 才会被打上 K48 泛素链交由蛋白酶体进行终结性降解12。此种时空分布上的 PTM 更替，完美地体现了免疫系统对 cGAS “该用时用、用完即毁”的精准治理逻辑。

### **7.2 乙酰化修饰对关键赖氨酸侧链的电荷剥夺与活性抑制**

除大分子修饰外，小分子的共价偶联同样能深刻改变泛素化图谱。赖氨酸在生理 pH 下带有正电荷，这不仅构成了 cGAS 捕获带负电荷 DNA 骨架的基础物理条件，也是其形成活性聚合物网络的前提。研究发现，组蛋白乙酰转移酶（如 KAT5）及其他多种调节酶类可介导 cGAS 的多位点乙酰化。最为值得关注的是，非甾体抗炎药物阿司匹林（Aspirin）被证实能够穿透细胞，直接乙酰化 cGAS 的 **K384、K394 和 K414** 位点25。

这三个位点的乙酰化带来了毁灭性的信号阻断：首先，乙酰基团的引入中和了赖氨酸原有的正电荷，使 cGAS 的结合界面不再亲和 DNA，无法形成缩合物；其次，回溯前文分析，K384 本是 RNF185 实施 K27 激活型泛素化的关键底物位点，而 K414 则是 p62 自噬降解与 USP14 解救争夺的枢纽位点22。乙酰化反应优先占据了这些残基的侧链氨基，从生化层面直接剥夺了后续泛素化级联发生的可能反应底物。这就意味着，乙酰化不仅作为独立的负性调节机制发挥作用，更通过直接“截胡”关键反应靶点，从根本上瘫痪了整个基于泛素化拓扑结构的 cGAS 活化与调控网络，强有力地阻断了过度自身免疫的病理恶化22。

## **8\. 基于 cGAS 泛素化图谱的靶向干预策略与转化医学前景**

将 cGAS 的泛素化网络从分子机理解析过渡到临床前转化医学，预示着一大批高特异性、低脱靶毒性的创新药物的诞生。

在针对**自身免疫与自身炎症性疾病**（如 SLE、AGS）的治疗开发中，除了传统的直接靶向 cGAS 催化口袋占据 ATP/GTP 位点的小分子抑制剂外，通过阻断上游激活型泛素化反应将是更具特异性的补充策略。开发小分子化合物变构抑制 TRIM56（针对 K347 单泛素化）或 RNF185（针对 K137/K384 的 K27 链），可从根源上平息异常的 cGAS 二聚化和过度敏感性12。另一方面，针对泛素化系统设计的药物也可以利用“刹车”机制，例如通过药物或代谢物间接增强 MARCH8 的表达与活性，加速对 cGAS K411 位点的空间位阻封锁，使天然免疫系统迅速降温41。

在**肿瘤免疫疗法与靶向放化疗**领域，cGAS 的泛素化生物学打开了两条截然不同却又殊途同归的干预路径： 其一，对于许多以免疫抑制为特征的“冷肿瘤”（Cold Tumors），肿瘤细胞极力维持极低的 cGAS 通路活性以规避 CD8+ T 细胞的免疫杀伤。此时，若能研发针对 DUBs 的药理学抑制剂——例如抑制 JOSD2，阻止其清除激活型的 K27 泛素链，从而阻止巨噬细胞向促肿瘤的 M2 表型极化56；或者针对降解途径，开发蛋白酶体抑制剂来延长细胞质中具备活性的 cGAS 的半衰期，即可在肿瘤微环境内诱导出高水平的 cGAMP，与现有的抗 PD-1/PD-L1 免疫检查点阻断疗法（ICB）形成强劲的协同抗肿瘤效应8。 其二，从**核内 cGAS 损害 DNA 修复机制**这一新视角出发，核内滞留的 cGAS 严重阻碍了癌细胞自身在发生损伤时的同源重组（HR）修复，从而改变了肿瘤的基因组进化方向并影响了其对放化疗的敏感度48。此时，针对 CRL5-SPSB3 降解复合体的药理调控成为前沿策略3。科学家正试图通过人工小分子底物招募技术（如蛋白降解靶向嵌合体，PROTAC），特异性利用或强化靶向 K427/K428 的泛素化机器，精准清除肿瘤细胞内的病理性 cGAS，从而重塑肿瘤基因组并打破化疗药物和 PARP 抑制剂遇到的疗效耐受瓶颈3。

## **9\. 结论与人源序列核心位点总结**

人类 cGAS 蛋白作为天然免疫监视的核心前哨，并非处于僵化的“开/关”状态，而是完全受控于一张由泛素-蛋白酶体系统、多维 E3 泛素连接酶网络及去泛素化酶（DUBs）交叉编织的动态翻译后修饰图谱。其发挥功能的确切后果，不仅取决于泛素链的拓扑形式（如单泛素、K27、K63、K48 和 K11），更被其表面特定的**赖氨酸（Lysine）残基**物理位置所严格界定。

基于对现存分子结构学、蛋白质组学及转化医学文献的全面梳理并消除物种间序列差异（人源与鼠源序列偏移），人源 cGAS 的特异性泛素化全景如下表所示：

| 人源 cGAS 赖氨酸位点 (Human cGAS) | 修饰酶 (E3 Ligase / DUB) | 泛素化连接类型 (Linkage Type) | 生化机制与生物学效应 (Molecular Mechanism & Biological Consequence) |
| :---- | :---- | :---- | :---- |
| **K137** | RNF185 (E3), JOSD2 (DUB) | K27-linked polyUb | K137 为非典型信号骨架起始点。多聚泛素化显著增强 cGAS 酶促活性（合成 cGAMP），促发强烈的抗病毒和自身免疫反应。JOSD2 可拮抗此作用。 |
| **K271** | USP29 (DUB) | K48-linked polyUb (去除) | 去除 K48 降解信号。组成性抑制 cGAS 蛋白在非感染期间的背景降解，维持受体在细胞质内的基本丰度阈值与免疫响应能力。 |
| **K347** (对应鼠 mK335) | TRIM56 (E3) | Mono-Ub | 由 TRIM56 特异性催化（Seo et al., 2018）。单泛素化发生于二聚化界面外围，提供构象推力，显著降低形成 cGAS-DNA 异源四聚体的能量壁垒，急剧放大底层免疫敏感性。**TRIM41 的 cGAS 泛素化位点尚未被实验鉴定。** |
| **K384** | RNF185 (E3), JOSD2 (DUB) | K27-linked polyUb | 与 K137 协同，作为非降解型信号脚手架，进一步稳定活性构象，促进 cGAMP 的高效生成。该位点同样为阿司匹林乙酰化的竞争靶标。 |
| **K411** | MARCH8 (E3) | K63-linked polyUb | 抑制性修饰。长链 K63 泛素紧邻 DNA 结合电荷面，产生巨大空间位阻，有效阻断 DNA 结合并关闭 cGAMP 生成，防止自身炎症发生。 |
| **K414** | 未知 E3, USP14 / TRIM14 (DUB) | K48-linked polyUb | K48 泛素链引导 p62 受体识别并启动溶酶体自噬降解程序；USP14 通过切割该链将 cGAS 从降解边缘挽救，延长其半衰期。同样受乙酰化竞争抑制。 |
| **K427 / K428** | CRL5-SPSB3 复合体 (E3) | K48-linked polyUb | 特异性针对核小体锚定的核内 cGAS。识别相邻 NN 基序以清除核内过剩受体，保障正常细胞周期运行和 DNA 同源重组的稳定性，防止慢性炎症。 |

综上所述，cGAS 的赖氨酸修饰网络是一部记录着机体病原体冲突与自身耐受平衡状态的高维“密码本”。通过解析并验证上述关键的赖氨酸残基，基础研究已经完成了从机制还原到干预靶点确证的历史跨越，为临床药理学界开辟靶向自身免疫病、神经炎症以及肿瘤微环境重塑的全新分子工具指明了明确的路径。

#### **Works cited**

1. Structural and Functional Analyses of DNA-Sensing and Immune Activation by Human cGAS, [https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0076983](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0076983)  
2. Structural and functional analyses of DNA-sensing and immune activation by human cGAS. | Literature citations | UniProt, [https://www.uniprot.org/citations/24116191](https://www.uniprot.org/citations/24116191)  
3. cGAS–STING: From immunology and oncology view \- PMC, [https://pmc.ncbi.nlm.nih.gov/articles/PMC12700745/](https://pmc.ncbi.nlm.nih.gov/articles/PMC12700745/)  
4. cGAS in nucleus: The link between immune response and DNA damage repair \- Frontiers, [https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2022.1076784/full](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2022.1076784/full)  
5. Roles of Emerging RNA-Binding Activity of cGAS in Innate Antiviral Response \- Frontiers, [https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2021.741599/full](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2021.741599/full)  
6. Post-Translational Modifications of cGAS-STING: A Critical Switch for Immune Regulation, [https://www.mdpi.com/2073-4409/11/19/3043](https://www.mdpi.com/2073-4409/11/19/3043)  
7. Research Advances in How the cGAS-STING Pathway Controls the Cellular Inflammatory Response \- Frontiers, [https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2020.00615/full](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2020.00615/full)  
8. Breaking barriers: The cGAS‐STING pathway as a novel frontier in cancer immunotherapy, [https://pmc.ncbi.nlm.nih.gov/articles/PMC12629869/](https://pmc.ncbi.nlm.nih.gov/articles/PMC12629869/)  
9. Lysine-11 ubiquitination drives type-I/III interferon induction by cGAS-STING and Toll-like receptors 3 and 4 \- PubMed, [https://pubmed.ncbi.nlm.nih.gov/41792265/](https://pubmed.ncbi.nlm.nih.gov/41792265/)  
10. The role of cGAS-STING pathway ubiquitination in innate immunity and multiple diseases, [https://pmc.ncbi.nlm.nih.gov/articles/PMC11868049/](https://pmc.ncbi.nlm.nih.gov/articles/PMC11868049/)  
11. Regulation and Function of the cGAS-STING Pathway: Mechanisms, Post-Translational Modifications, and Therapeutic Potential in Immunotherapy \- PMC, [https://pmc.ncbi.nlm.nih.gov/articles/PMC11911240/](https://pmc.ncbi.nlm.nih.gov/articles/PMC11911240/)  
12. The Regulation of cGAS \- PMC, [https://pmc.ncbi.nlm.nih.gov/articles/PMC5934468/](https://pmc.ncbi.nlm.nih.gov/articles/PMC5934468/)  
13. Regulation and inhibition of the DNA sensor cGAS \- PMC \- NIH, [https://pmc.ncbi.nlm.nih.gov/articles/PMC7726805/](https://pmc.ncbi.nlm.nih.gov/articles/PMC7726805/)  
14. Regulation and Function of the cGAS-STING Pathway: Mechanisms, Post-Translational Modifications, and Therapeutic Potential in Immunotherapy \- PubMed, [https://pubmed.ncbi.nlm.nih.gov/40098909/](https://pubmed.ncbi.nlm.nih.gov/40098909/)  
15. cGAS-Stimulator of Interferon Genes Signaling in Central Nervous System Disorders \- Aging and disease, [https://www.aginganddisease.org/EN/article/downloadArticleFile.do?attachType=PDF\&id=148070](https://www.aginganddisease.org/EN/article/downloadArticleFile.do?attachType=PDF&id=148070)  
16. Primary disorders of polyubiquitination: Dual roles in autoinflammation and immunodeficiency \- Rockefeller University Press, [https://rupress.org/jem/article/222/5/e20241047/277388/Primary-disorders-of-polyubiquitination-Dual-roles](https://rupress.org/jem/article/222/5/e20241047/277388/Primary-disorders-of-polyubiquitination-Dual-roles)  
17. E3 ubiquitin ligases: styles, structures and functions \- PMC, [https://pmc.ncbi.nlm.nih.gov/articles/PMC8607428/](https://pmc.ncbi.nlm.nih.gov/articles/PMC8607428/)  
18. Ubiquitin \- Wikipedia, [https://en.wikipedia.org/wiki/Ubiquitin](https://en.wikipedia.org/wiki/Ubiquitin)  
19. Survey of the human proteostasis network: the ubiquitin-proteasome system | bioRxiv, [https://www.biorxiv.org/content/10.64898/2026.03.13.711689v1.full](https://www.biorxiv.org/content/10.64898/2026.03.13.711689v1.full)  
20. Deubiquitinase OTULIN dampens RIG-I-dependent antiviral signaling by removing linear ubiquitination from TRAF6 | PNAS, [https://www.pnas.org/doi/10.1073/pnas.2517201123](https://www.pnas.org/doi/10.1073/pnas.2517201123)  
21. The role of cGAS-STING pathway ubiquitination in innate immunity and multiple diseases, [https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2025.1522200/full](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2025.1522200/full)  
22. Post-Translational Modifications of Proteins in Cytosolic Nucleic Acid Sensing Signaling Pathways \- PMC, [https://pmc.ncbi.nlm.nih.gov/articles/PMC9250978/](https://pmc.ncbi.nlm.nih.gov/articles/PMC9250978/)  
23. Research Advances in How the cGAS-STING Pathway Controls the Cellular Inflammatory Response \- PMC, [https://pmc.ncbi.nlm.nih.gov/articles/PMC7198750/](https://pmc.ncbi.nlm.nih.gov/articles/PMC7198750/)  
24. Human cGAS-DNA Complex Structure Insights | PDF \- Scribd, [https://www.scribd.com/document/859335645/2018-Cell-Zhou-Et-Al-Structure-of-the-Human-CGAS-DNA-Complex-Reveals](https://www.scribd.com/document/859335645/2018-Cell-Zhou-Et-Al-Structure-of-the-Human-CGAS-DNA-Complex-Reveals)  
25. cGAS–STING, an important signaling pathway in diseases and their therapy \- PMC \- NIH, [https://pmc.ncbi.nlm.nih.gov/articles/PMC10960729/](https://pmc.ncbi.nlm.nih.gov/articles/PMC10960729/)  
26. E3 ubiquitin-protein ligase TRIM56 \- Homo sapiens (Human) | UniProtKB | UniProt, [https://www.uniprot.org/uniprotkb/Q9BRZ2/entry](https://www.uniprot.org/uniprotkb/Q9BRZ2/entry)  
27. The DNA Sensor cGAS is Decorated by Acetylation and Phosphorylation Modifications in the Context of Immune Signaling \- PMC, [https://pmc.ncbi.nlm.nih.gov/articles/PMC7338091/](https://pmc.ncbi.nlm.nih.gov/articles/PMC7338091/)  
28. CGAS \- Cyclic GMP-AMP synthase \- Homo sapiens (Human) | UniProtKB | UniProt, [https://www.uniprot.org/uniprotkb/Q8N884/entry](https://www.uniprot.org/uniprotkb/Q8N884/entry)  
29. E3 ubiquitin-protein ligase TRIM41 \- Homo sapiens (Human) | UniProtKB | UniProt, [https://www.uniprot.org/uniprotkb/Q8WV44/entry](https://www.uniprot.org/uniprotkb/Q8WV44/entry)  
30. TRIM41 Gene \- GeneCards | TRI41 Protein, [https://www.genecards.org/cgi-bin/carddisp.pl?gene=TRIM41](https://www.genecards.org/cgi-bin/carddisp.pl?gene=TRIM41)  
31. RINCK-mediated monoubiquitination of cGAS promotes antiviral innate immune responses, [https://pmc.ncbi.nlm.nih.gov/articles/PMC5944131/](https://pmc.ncbi.nlm.nih.gov/articles/PMC5944131/)  
32. Nuclear cGAS restricts L1 retrotransposition by promoting TRIM41-mediated ORF2p ubiquitination and degradation \- ResearchGate, [https://www.researchgate.net/publication/376453578\_Nuclear\_cGAS\_restricts\_L1\_retrotransposition\_by\_promoting\_TRIM41-mediated\_ORF2p\_ubiquitination\_and\_degradation](https://www.researchgate.net/publication/376453578_Nuclear_cGAS_restricts_L1_retrotransposition_by_promoting_TRIM41-mediated_ORF2p_ubiquitination_and_degradation)  
33. cGAS- Stimulator of Interferon Genes Signaling in Central Nervous System Disorders \- PMC, [https://pmc.ncbi.nlm.nih.gov/articles/PMC8460300/](https://pmc.ncbi.nlm.nih.gov/articles/PMC8460300/)  
34. K11-Linked Polyubiquitination in Cell Cycle Control Revealed by a K11 Linkage-Specific Antibody | Request PDF \- ResearchGate, [https://www.researchgate.net/publication/45286197\_K11-Linked\_Polyubiquitination\_in\_Cell\_Cycle\_Control\_Revealed\_by\_a\_K11\_Linkage-Specific\_Antibody](https://www.researchgate.net/publication/45286197_K11-Linked_Polyubiquitination_in_Cell_Cycle_Control_Revealed_by_a_K11_Linkage-Specific_Antibody)  
35. The ubiquitin–UBD network.Ubiquitin can be covalently attached to... \- ResearchGate, [https://www.researchgate.net/figure/The-ubiquitin-UBD-networkUbiquitin-can-be-covalently-attached-to-target-proteins-as-a\_fig2\_26831034](https://www.researchgate.net/figure/The-ubiquitin-UBD-networkUbiquitin-can-be-covalently-attached-to-target-proteins-as-a_fig2_26831034)  
36. Scientists identify ANKIB1 as key regulator of innate immune signaling \- News-Medical.Net, [https://www.news-medical.net/news/20260306/Scientists-identify-ANKIB1-as-key-regulator-of-innate-immune-signaling.aspx](https://www.news-medical.net/news/20260306/Scientists-identify-ANKIB1-as-key-regulator-of-innate-immune-signaling.aspx)  
37. Lysine-11 ubiquitination drives type-I/III interferon induction by cGAS-STING and Toll-like receptors 3 and 4\. \- Apollo, [https://www.repository.cam.ac.uk/items/ae90bcda-27d1-4f9f-8955-19809523507d](https://www.repository.cam.ac.uk/items/ae90bcda-27d1-4f9f-8955-19809523507d)  
38. MARCH8-mediated ubiquitination regulates expression of the antiviral protein IFITM3 \- PMC, [https://pmc.ncbi.nlm.nih.gov/articles/PMC12702058/](https://pmc.ncbi.nlm.nih.gov/articles/PMC12702058/)  
39. MARCHF8 Gene \- GeneCards | MARH8 Protein, [https://www.genecards.org/cgi-bin/carddisp.pl?gene=MARCHF8](https://www.genecards.org/cgi-bin/carddisp.pl?gene=MARCHF8)  
40. The Membrane-Associated MARCH E3 Ligase Family: Emerging Roles in Immune Regulation \- Frontiers, [https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2019.01751/full](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2019.01751/full)  
41. MARCH8 attenuates cGAS-mediated innate immune responses through ubiquitylation \- PubMed, [https://pubmed.ncbi.nlm.nih.gov/35503863/](https://pubmed.ncbi.nlm.nih.gov/35503863/)  
42. cGAS, an innate dsDNA sensor with multifaceted functions \- PMC, [https://pmc.ncbi.nlm.nih.gov/articles/PMC12135389/](https://pmc.ncbi.nlm.nih.gov/articles/PMC12135389/)  
43. (PDF) MARCH8 attenuates cGAS-mediated innate immune responses through ubiquitylation \- ResearchGate, [https://www.researchgate.net/publication/360354118\_MARCH8\_attenuates\_cGAS-mediated\_innate\_immune\_responses\_through\_ubiquitylation](https://www.researchgate.net/publication/360354118_MARCH8_attenuates_cGAS-mediated_innate_immune_responses_through_ubiquitylation)  
44. MARCHF8 profile page \- Open Targets Platform, [https://platform.opentargets.org/target/ENSG00000165406](https://platform.opentargets.org/target/ENSG00000165406)  
45. E3 ubiquitin-protein ligase MARCHF8 \- Homo sapiens (Human) | UniProtKB | UniProt, [https://www.uniprot.org/uniprotkb/Q5T0T0/entry](https://www.uniprot.org/uniprotkb/Q5T0T0/entry)  
46. USP29 maintains the stability of cGAS and promotes cellular antiviral responses and autoimmunity \- PMC, [https://pmc.ncbi.nlm.nih.gov/articles/PMC7608407/](https://pmc.ncbi.nlm.nih.gov/articles/PMC7608407/)  
47. TRIM14 Inhibits cGAS Degradation Mediated by Selective Autophagy Receptor p62 to Promote Innate Immune Responses | Request PDF \- ResearchGate, [https://www.researchgate.net/publication/308571988\_TRIM14\_Inhibits\_cGAS\_Degradation\_Mediated\_by\_Selective\_Autophagy\_Receptor\_p62\_to\_Promote\_Innate\_Immune\_Responses](https://www.researchgate.net/publication/308571988_TRIM14_Inhibits_cGAS_Degradation_Mediated_by_Selective_Autophagy_Receptor_p62_to_Promote_Innate_Immune_Responses)  
48. SPSB3 targets nuclear cGAS a, Results from the siRNA screen... \- ResearchGate, [https://www.researchgate.net/figure/SPSB3-targets-nuclear-cGAS-a-Results-from-the-siRNA-screen-highlighting-cGAS-GFP-nuclear\_fig3\_378550542](https://www.researchgate.net/figure/SPSB3-targets-nuclear-cGAS-a-Results-from-the-siRNA-screen-highlighting-cGAS-GFP-nuclear_fig3_378550542)  
49. The Odyssey of cGAS: from Cytosol to Nucleus \- PMC, [https://pmc.ncbi.nlm.nih.gov/articles/PMC11542052/](https://pmc.ncbi.nlm.nih.gov/articles/PMC11542052/)  
50. 1 A p97-PML NBs axis is essential for chromatin-bound cGAS untethering and degradation upon senescence-prone DNA damage Florent \- bioRxiv, [https://www.biorxiv.org/content/10.1101/2025.10.24.684316.full.pdf](https://www.biorxiv.org/content/10.1101/2025.10.24.684316.full.pdf)  
51. The CRL5–SPSB3 ubiquitin ligase targets nuclear cGAS for degradation \- ResearchGate, [https://www.researchgate.net/publication/378550542\_The\_CRL5-SPSB3\_ubiquitin\_ligase\_targets\_nuclear\_cGAS\_for\_degradation](https://www.researchgate.net/publication/378550542_The_CRL5-SPSB3_ubiquitin_ligase_targets_nuclear_cGAS_for_degradation)  
52. Targeting the PRMT1-cGAS-STING signaling pathway to enhance the anti-tumor therapeutic efficacy \- PMC, [https://pmc.ncbi.nlm.nih.gov/articles/PMC11867627/](https://pmc.ncbi.nlm.nih.gov/articles/PMC11867627/)  
53. USP29 maintains the stability of cGAS and promotes cellular antiviral responses and autoimmunity | Request PDF \- ResearchGate, [https://www.researchgate.net/publication/341649766\_USP29\_maintains\_the\_stability\_of\_cGAS\_and\_promotes\_cellular\_antiviral\_responses\_and\_autoimmunity](https://www.researchgate.net/publication/341649766_USP29_maintains_the_stability_of_cGAS_and_promotes_cellular_antiviral_responses_and_autoimmunity)  
54. USP27X \- Ubiquitin carboxyl-terminal hydrolase 27 \- Homo sapiens (Human) | UniProtKB, [https://www.uniprot.org/uniprotkb/A6NNY8/entry](https://www.uniprot.org/uniprotkb/A6NNY8/entry)  
55. 389856 \- Gene ResultUSP27X ubiquitin specific peptidase 27 X-linked \[ (human)\] \- NCBI, [https://www.ncbi.nlm.nih.gov/gene/389856](https://www.ncbi.nlm.nih.gov/gene/389856)  
56. Full article: Deubiquitinating enzyme JOSD2 modulates cGAS to facilitate immune evasion in colorectal cancer \- Taylor & Francis, [https://www.tandfonline.com/doi/full/10.1080/2162402X.2025.2590245](https://www.tandfonline.com/doi/full/10.1080/2162402X.2025.2590245)  
57. (PDF) Deubiquitinating enzyme JOSD2 modulates cGAS to facilitate immune evasion in colorectal cancer \- ResearchGate, [https://www.researchgate.net/publication/398409656\_Deubiquitinating\_enzyme\_JOSD2\_modulates\_cGAS\_to\_facilitate\_immune\_evasion\_in\_colorectal\_cancer](https://www.researchgate.net/publication/398409656_Deubiquitinating_enzyme_JOSD2_modulates_cGAS_to_facilitate_immune_evasion_in_colorectal_cancer)