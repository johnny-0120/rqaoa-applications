# Recursive QAOA 應用方向研究

> 2026-11-12後續實測：這份文件保留測試前的探索推論。已執行的點雲比較未見RQAOA優勢，改善主要來自幾何模型；最新證據見 [實測報告](<./output/點雲法向量_RQAOA與古典方法比較_2026-11-12.md>)。

研究跨日：2026-11-11 至 11-12（Asia/Taipei）。作者：Codex AI。這是文獻查核與應用提案，尚非 Johnny 採用的選題，未執行本報告提出的應用實驗。

## 研究結論

優先考慮 **3D 點雲局部面片的法向量方向校正**，其次是 **具雜訊及缺測的感測器網路極性校正**。兩者的未知量本來就是正／反向，且可靠區域可以逐步合併，因此與 RQAOA 的「估計關係 → 固定相對符號 → 消去變數 → 解剩餘核心 → 還原」機制契合。

本次公開網路檢索未找到這兩項特定用途的直接 RQAOA 研究；這是初步研究缺口，不能保證全球無人做過。二者的問題建模、圖同步、古典縮圖與領域求解方法已有大量研究。可爭取的是具體應用、可信關係的定義及端到端實驗貢獻，不能把同一個 Ising 模型換名字就稱為新演算法。

另外保留兩個較次要方向：分散式二分標籤對齊、基因定相區塊拼接。後者已有很近的量子啟發式 MaxCut 方法，選題風險較高。

## 從哪些特色推導應用

| RQAOA 特色 | 適合找的問題結構 | 必須保留的限制 |
|---|---|---|
| 使用二點期望值 Mij=〈ZiZj〉，可固定同向或反向關係 | 相對方向比絕對標籤容易知道；正反號本來就是物理或資料意義 | Mij 是 QAOA 狀態的相關性，不是資料相關係數，也不是「該關係正確」的後驗機率 |
| 固定 zj=s zi 後，精確更新剩餘目標 | 一批變數可被代表變數吸收，局部區塊具有共同符號 | 代入在所選限制下是精確的；選哪個限制是啟發式，可能排除全域最佳解 |
| 透過反覆縮圖，跨步整合非局部資訊 | 局部判斷容易被矛盾迴圈誤導，需要整體一致性 | 非局部處理也有古典方法；不等於量子獨有能力 |
| 最後解小型剩餘問題並回填 | 大部分區域容易，少數接縫或橋接邊困難 | 小核心求得最佳，不代表原問題全域最佳 |
| 各輪可維持小的 QAOA 層數 p | 願意以多輪處理換取較少量子演化層數 | 起始仍通常需要 n 個邏輯 qubits；p 固定不代表編譯電路深度固定，縮圖可使圖變稠密 |
| 保留消去關係及還原紀錄 | 需要追查哪一步造成錯誤或如何形成群組 | 縮圖樹只能說明演算法過程，不能當正確性或物理因果證明 |

以上機制以 Bravyi 等原始論文附錄 C 與 OpenQAOA 官方文件為依據。[S1][S2] 應用適配與下面的排序是本次 AI 推論。

共同模型可寫成：

\[
\min_{z_i\in\{-1,+1\}} L(z)=\sum_{(i,j)} w_{ij}\frac{1-r_{ij}z_i z_j}{2},\quad w_{ij}\geq0,\ r_{ij}\in\{-1,+1\}.
\]

r 表示資料提供的相對符號；w 表示資料證據強度。這和演算法內部的 Mij 是不同量。若以 zj=s zi 消去 j，在一般 Ising 能量 E=c+Σhi zi+ΣJij zi zj 下，更新為 c'=c+sJij、h'i=hi+s hj、J'ik=Jik+sJjk，其餘不變。更新後仍為二次模型。

沒有絕對參考時，z 與 -z 等價。無外場、標準對稱初態與 mixer 下，〈Zi〉也可能恆為零，因此應利用二點關係而非依賴單點符號。每個不連通分量都各有一個無法辨識的整體翻轉自由度。固定參考只是在選座標；要判斷真正向外或正極，仍需外部證據。

## 方向一：3D 點雲面片的法向量方向校正（首選）

**問題。** 3D 掃描可估計每點的法線軸，卻常不知道箭頭應指向哪邊。若局部面片已內部一致，只剩每片整體保留或翻轉的選擇。錯誤方向會影響後續曲面重建。局部／全局分階段和面片方向傳播已有研究，並非本次發明。[S3][S4]

**從特色發想。** 每片一個 z；RQAOA 固定的 zi=zj 或 zi=-zj，恰好可解讀為「兩片維持相對方向」或「其中一片翻轉」。先吸收容易判斷的面片，再處理狹縫、遮擋、薄殼附近互相矛盾的接縫。這比把每個原始點都硬塞進電路更適合作為小型研究，但面片化是額外古典前處理，不是 RQAOA 自動突破 qubit 限制。

**模型。** 從面片重疊區／可信鄰接取得 r 與 w，套用共同模型。對平滑且確認在同一側曲面上的相鄰樣本，可用法線內積符號作起點；不能在薄板兩側或尖角直接用「空間很近就同向」，否則求解器再強也會解錯建模。可利用掃描視點、可見性或幾何分隔建立較可靠邊，另保留不使用這些資訊的消融對照。

**最小實驗。** 先以有真值的球、帶洞形狀、薄殼／雙層近接面生成點雲，分成 12、16、20 個局部一致面片，隨機翻轉整片並加入缺邊、離群點及錯誤鄰接。再進入公開真實掃描；資料授權與現成程式尚未核驗，不能說已下載。所有方法必須用同一份面片及邊，另做端到端管線比較。

**必要對照。** 最大可靠性生成樹傳播、譜同步、SDP、模擬退火／多起點局部搜尋、同樣縮圖但改用古典關係分數的版本。端到端另對照 dipole propagation、WNNC 或可取得的現代法向定向方法，不能只挑弱 MST。[S3][S4][S5]

**指標。** 去除整體翻轉歧義後的面片／點定向錯誤率、面積加權錯誤、原始圖代價、曲面重建品質及總時間。能量較低不保證真值方向較準；必須同時回報。

**可研究的具體缺口。** 幾何可靠性分級後，RQAOA 的消去次序是否能減少「錯誤跨薄殼連邊」的污染？保留弱幾何邊至最後處理是否有效？這些是待驗證的應用策略；可靠度加權、延後決策或回溯本身不能聲稱原創。

**淘汰條件。** 同預算下始終不如譜同步／SDP／領域方法；收益全來自更好的前處理；或簡單局部方向目標無法對應曲面品質。此時保留診斷結果，不包裝成量子優勢。

**新穎性狀態。** 未找到直接 RQAOA × 面片法向量定向研究；已有成熟的古典全局定向與面片傳播研究。此方向最適合展示「RQAOA 機制與問題變數一致」，尚不能保證性能。

## 方向二：感測器網路的極性校正（次選）

**問題。** 某些通道因接線、安裝或處理流程而整體反號。兩個通道可以估計相對正負，但不一定知道誰才是反接者。EarthScope 的實際品質指標就利用相鄰測站交叉相關檢查相對極性，並明示單一站對不能定位哪一站反向。[S6]

**從特色發想。** 每個通道一個 z；原始資料有很多局部相對證據，卻可能受雜訊、缺測與延遲誤差影響，形成矛盾迴圈。RQAOA 先建立通道間的相對群組，最後透過已知參考通道還原絕對符號。這裡解的是設備通道的二元極性，不是任意角度、時間同步或整個連續校準問題。

**模型與假設。** 在具有共同校準訊號、已處理時間延遲及已知幾何差異的條件下，利用多個事件／時間窗估計 r、w，再解共同模型。不能把所有負相關都當反接：真正物理響應也可能反相。應先用人工注入的極性錯誤建立可驗證真值，實際場域再審核物理條件。

**最小實驗。** 12–20 通道，生成共同訊號、不同增益／延遲與雜訊，注入已知反接，再掃描觀測圖連通性、錯邊比例與 SNR。資料估計 r 的訓練窗和驗證訊號必須分開。先做受控振動或音訊校準；地震資料只作領域延伸，不能直接移植簡化物理假設。

**必要對照與指標。** 生成樹、譜同步、SDP、局部搜尋／退火及古典縮圖；回報極性恢復率、已知參考下的絕對錯誤率、未見事件的相干疊加增益、總處理時間。圖同步本身已有強古典方法，須納入。[S5]

**可研究的具體缺口。** 當可靠觀測集中在若干群組，群組間只有少量互相衝突的橋接邊時，RQAOA 是否比一次譜取整更能恢復相對極性？若不行，縮圖紀錄是否能指出跨群組錯誤如何擴散？前者是效能假說，後者是診斷用途，都未經本輪實驗。

**淘汰條件。** 物理反相與接線反相不可辨識；圖不連通又沒有每分量參考；或強古典方法同成本已全部解對。不能為製造困難而只保留人造極端噪音，再宣稱實際適用。

**新穎性狀態。** 未找到直接 RQAOA × 設備通道極性校正研究；相對極性估計及同步模型不是新問題。

## 方向三：分散式二分標籤／局部分割結果對齊（備選）

**問題。** 多個區域／節點各做二分群，A/B 的名字可以任意交換。即使每份局部結果不錯，重疊處仍可能因標籤反轉而不一致。例如各影像區塊各自做「兩種材質」分群，再對齊共享像素。這是沒有固定類別語意的二分輸出；若已有固定「前景／背景」語意，不應憑空製造可交換性。

**適配與模型。** 每份局部結果一個 z，重疊樣本提供兩份結果應同名或換名的投票 r 與權重 w。RQAOA 合併的是整份局部結果的標籤座標，不是樣本本身，也不是讓量子電腦訓練模型。

**最小實驗。** 12–20 個重疊子資料集，控制重疊比例、局部分群錯誤與圖上的矛盾迴圈；測對齊後的共享樣本一致性，以及在具真值資料上的最終分群品質。比較各自對齊單一參考、生成樹、譜同步、SDP 與古典縮圖。

**候選貢獻與風險。** 研究「重疊圖的弱橋接」是否造成不可逆錯誤，及 RQAOA 是否有可重現改善。這本質上仍是二元同步；與前兩方向共享數學，不能算三種新演算法。多於兩類會變成置換同步，不能直接沿用一個 Ising bit 的模型。重疊無資訊時無法辨識。

**新穎性狀態。** 未找到完全匹配「RQAOA × 局部二分結果座標對齊」的直接研究；尚未完整審查 ensemble clustering、permutation synchronization、divide-and-conquer QAOA 的全部交叉文獻，因此新穎性信心低於前兩項。暫作備選，不推薦直接宣稱首創。

## 方向四：基因定相區塊的正反向拼接（高重疊風險備選）

**問題。** 已完成的局部二倍體定相區塊仍可能整體交換兩套單倍型標籤。跨區塊 reads 或其他長距離證據可提供相對一致／交換關係；候選是一個區塊一個 z，而非一開始就做整條染色體每個變異位點的量子求解。

**適配。** 合併兩個區塊的正／反向標籤，正好對應相對符號消去，且長距離弱證據形成困難核心。但跨區塊來源及權重需符合基因資料模型，簡單投票模型不等於完整 minimum-error-correction。

**已知重疊。** QHap 已將定相寫成 MaxCut，使用 GPU 量子退火啟發式的 simulated bifurcation；因此「量子概念＋定相＋MaxCut」不能稱新。[S7] 本次未找到直接 RQAOA 區塊拼接研究，不代表只是替換求解器就有足夠貢獻。

**最小實驗與對照。** 先合成 12–20 個具有已知方向的區塊，加入 bridge read 錯誤與缺測；比較加權 MST、譜同步、SDP、局部搜尋，後續必須對照 QHap 及成熟領域工具。回報 switch error、正確拼接區塊比例、連續度與時間，不能只看 MaxCut 分數。工具部署與資料授權本輪未驗證。

**判斷。** 生物資料處理與強對照成本較高，且已存在近鄰研究，不列競賽首選。只有掌握具體跨區塊困難與領域資料時才值得投入。

## 已查到前人工作，應排除的泛用「創新」說法

| 想法 | 直接或近鄰先例 | 本次處理 |
|---|---|---|
| RQAOA 用於汽車噴漆排程 | Optimization-Free RQAOA for Binary Paint Shop [S8] | 已有直接應用 |
| RQAOA 用於無線頻道／干擾資源配置 | Chen 等 [S9] | 已有直接應用 |
| RQAOA 加 RL 選消去關係 | RL-RQAOA [S10] | 已有人做，不能當核心新意 |
| 看困難度分配 shots | Adaptive Shot Allocation [S11] | 已有直接研究 |
| 選錯回溯、保留符號分支 | QIRO [S12]、11/5 新預印本 beam RQAOA [S13] | 可作對照／工具，不能直接聲稱原創 |
| 直接換成投資組合題目 | beam RQAOA 已包含基數限制投資組合 [S13] | 本次不列新用途 |
| 只因 QUBO 可表達就推薦任何排程／物流 | 可表達不等於適配；硬限制／輔助變數會增加難度 | 先選相對符號本來就有意義的問題 |

這份表不是完整應用綜述，也沒有把未全文核讀的摘要成果當已獨立重現。

## 建議如何銜接既有實作

本輪讀到本機 `output/第9題_RQAOA研究回報_團隊版.md`（文件日期11/7），其自述208案例中，RQAOA找到最佳解197例、古典局部搜尋206例、SDP/GW 207例，且成本未嚴格匹配。本次未重跑，不把這些數字當新驗收，也不歸為已證實的應用結果。它支持一個務實判斷：不能再只以「能解圖」當主要貢獻。

已查看 `rqaoa_experiment.py` 的 p=1 成對相關公式與 `rqaoa` 入口，具備加權圖基礎；但原註解明示參數搜尋採整數權重假設。新應用的實數、有正負權重及可能的外場，必須另外核對角度範圍、正負號、offset 與還原計分，不能宣稱現程式已支援全部新問題。原 MaxCut 最大化與本報告最小化能量的符號、常數項也需統一。

**建議的最小研究問題：**「面片幾何觀測含錯誤跨面連邊時，RQAOA 是否能在相同總資源下改善全局法向量定向？它在哪一步把局部幾何錯誤擴散為全域方向錯誤？」

若採首選，實驗先分兩層：固定同一份幾何圖比較求解器，再比較整個重建管線。否則無法分清收益是幾何建模、古典面片化或 RQAOA 帶來。

共同實驗規格建議如下（AI建議，未執行、未承諾期限）：

1. 先做可枚舉的小圖，保留未見形狀／未見雜訊條件作測試，避免用測試集反覆選規則。方法共用圖、前處理、真值與評分。
2. 固定深度、cutoff、shots、重啟次數，回報全部遞迴步驟的電路評估、總 shots、雙 qubit 閘及編譯深度、古典時間及端到端時間。有限shots應來自共同 bitstring 樣本，保留相關估計的聯合結構。
3. 同時比較標準 QAOA 的固定預算最好樣本、原 RQAOA、古典縮圖、譜／SDP、退火／局部搜尋。量子方法與古典方法分別掃預算，畫品質－成本曲線，避免把不同資源混成公平排名。
4. 小圖用枚舉測原始最佳值及每步限制造成的損失；大圖以最佳已知值與可用界限回報，不假裝知道真最佳。
5. 主指標用原始錯誤率／代價差；有正負權重或常數平移時，不直接用可能失真的「能量近似比」。需另報中位數、分布及跨實例不確定性。
6. 對無噪音可一致滿足的圖先做 sanity check，真正比較放在有領域意義的矛盾迴圈；不能用純樹或可直接傳播的情況宣稱量子優勢。

原始論文已給出 Ising p=1 相關量的有效古典計算方法。[S1] 因此先用 p=1 找應用結構很合理，但不構成量子計算優勢證據；增加 p 也不能自動保證難以模擬或勝過古典。更適合先主張可重現的混合求解／失敗機制研究，再依證據決定是否值得上真機。

## 檢索工作紀錄與證據界線

方法：先讀 Vault 的 Graphify INTEGRATION 與相關技能，從圖譜定位競賽專案，讀原始筆記；競賽原文 SHA256 與 sync-state 相符。圖譜可見詞彙含 QAOA、R-QAOA與Warm-start QAOA分線。再直接查原始論文、作者預印本及官方文件。一般查詢未重建全庫。

檢索涵蓋名稱 Recursive QAOA、RQAOA、recursive quantum approximate optimization，並加入應用詞及近鄰方法 QAOA、quantum、Ising、synchronization。主要實際查詢包含：

- `recursive QAOA applications correlation rounding optimization RQAOA`
- `"recursive QAOA" applications 2025 2026`
- `"recursive QAOA" "applications"`
- `"recursive QAOA" "normal" orientation`
- `"QAOA" "surface normal" orientation`
- `QAOA normal orientation point cloud Ising`
- `"QAOA" "normal orientation"`
- `"normal orientation" "quantum"`
- `"normal orientation" "Ising"`
- `"RQAOA" "point cloud"`
- `"QAOA" "normal" "point clouds"`
- `"recursive quantum approximate optimization" "normal"`
- `"QAOA" "sensor polarity"`
- `"QAOA" "polarity" calibration`
- `quantum optimization sensor polarity synchronization Ising`
- `"QAOA" "Z2" synchronization`
- `"recursive QAOA" "synchronization"`
- `"recursive quantum approximate optimization" "polarity"`
- `"quantum" "polarity correction" sensor`
- `microphone polarity calibration sign ambiguity synchronization paper`
- `seismic sensor polarity reversal detection cross correlation paper`
- `"QAOA" "label alignment"`
- `"QAOA" "label synchronization"`
- `"recursive QAOA" "label" alignment segmentation`
- `"RQAOA" "ensemble"`
- `"binary" "cluster ensemble" "sign" synchronization`
- `"RQAOA" "haplotype"`
- `recursive QAOA haplotype phasing genome`
- `"QAOA" "haplotype" phasing`
- `"recursive QAOA" "haplotype phasing"`
- `"recursive QAOA" "haplotype" phasing block`

部分精確查詢只回傳無關內容；這只能記為沒有取得直接命中，不能當不存在的證明。沒有全面逐庫查 Scopus／Web of Science、學位論文、專利、未公開投稿與所有代碼庫；新穎性屬初查，仍需追引用、換同義詞、查近期全文與專利。未取得結果的查詢不補造論文。

閱讀深度：S1已讀摘要及附錄C相關機制／古典可模擬性段落；S2、S6讀官方操作／方法正文；S3–S5、S7–S13以原始論文摘要／書目為主，足以辨識鄰近研究，不足以驗收其數據或完整排除方法重疊。文中未直接引用大段原文。

本輪完成的是選題研究與保存，沒有新應用資料集、應用效能結果或量子真機實驗。後續採不採用、團隊題目與官方參賽期限都未在本輪定案。

## 可回查來源

- [S1] Bravyi et al., Obstacles to State Preparation and Variational Optimization from Symmetry Protection，2019預印本／PRL 2020；RQAOA原始研究。[摘要](https://arxiv.org/abs/1910.08980) · [附錄C全文](https://arxiv.org/html/1910.08980v1)
- [S2] OpenQAOA 官方 Recursive QAOA workflow。[文件](https://openqaoa.entropicalabs.com/workflows/recursive-qaoa/)
- [S3] Metzer et al., Orienting Point Clouds with Dipole Propagation，SIGGRAPH 2021；局部面片與全局定向。[原始論文](https://arxiv.org/abs/2105.01604)
- [S4] Lin et al., Fast and Globally Consistent Normal Orientation based on the Winding Number Normal Consistency，2024；全局法向量強對照。[原始論文](https://arxiv.org/abs/2405.16634)
- [S5] Singer, Angular Synchronization by Eigenvectors and Semidefinite Programming，2009預印本／2011期刊；同步、譜及SDP方法。[原始論文](https://arxiv.org/abs/0905.3174)
- [S6] EarthScope MUSTANG polarity_check；相鄰通道交叉相關、相對極性及辨識限制。[官方文件](https://service.earthscope.org/mustang/metrics/docs/1/desc/polarity_check/)
- [S7] Zhang et al., QHap: Quantum-Inspired Haplotype Phasing，2026；MaxCut與GPU simulated bifurcation。[原始預印本](https://arxiv.org/abs/2603.25762)
- [S8] Mooney et al., Optimisation-Free Recursive QAOA for the Binary Paint Shop Problem，2025預印本、2026修訂；直接應用先例。[原始論文](https://arxiv.org/abs/2507.10908)
- [S9] Chen et al., Recursive QAOA for Interference-Aware Resource Allocation in Wireless Networks，2026。[原始預印本](https://arxiv.org/abs/2602.07483)
- [S10] Patel et al., Reinforcement Learning Assisted Recursive QAOA，2022預印本／2024期刊。[原始論文](https://arxiv.org/abs/2207.06294)
- [S11] Lee & Kim, Adaptive Shot Allocation for Recursive QAOA via Reinforcement Learning，2026。[原始預印本](https://arxiv.org/abs/2605.26544)
- [S12] Finžgar et al., Quantum-Informed Recursive Optimization Algorithms，2023預印本／2024期刊；問題特定縮減與回溯。[原始論文](https://arxiv.org/abs/2308.13607)
- [S13] Ray & Sharma, Sign-hedged Beam QAOA: Scaling Quantum Optimization Algorithms，2026-11-05；兩種符號分支及投資組合案例。[原始預印本](https://arxiv.org/abs/2611.08646)

## 本次需求原話

來源：Johnny，本次 Codex 對話；研究開始11/11、完成保存11/12。

> 請你幫我找找看 Recursive QAOA 有什麼可以應用的方向。請找那些還沒有人做過、但它的特性很適合做的方向，從它的特色去發想，看看應用可以做哪一些方面？並且把工作記錄存到 Obsidian 的量子競賽那一個專案裡面

本文件將研究建議與既有事實分開；沒有把建議轉成 Johnny 決定。
