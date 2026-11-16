"""Concise team-facing progress report; preserve completed/proposed distinction."""
from pathlib import Path
import json,re,html
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'output';PDF=OUT/'pdf'
A=json.loads((ROOT/'research/q9_results/analysis.json').read_text(encoding='utf-8'))
pages=[('第 9 題：Recursive QAOA 初步研究回報','''
回報人：Johnny｜日期：2026-11-07

## 目前判斷

我負責驗證 Recursive QAOA（R-QAOA）能否改善第 9 題中，標準 QAOA 在不同圖結構上的表現。目前程式與比較實驗已完成，但結果仍屬可行性驗證與失敗分析，尚不足以作為參賽的主要突破。

R-QAOA 在這批鏈、環、星狀、樹、完全圖及二分圖上都能找到最佳解；但社群圖、一般隨機圖與部分控制圖仍會失敗。加入重複抽樣及古典對照後，沒有看到它全面領先的證據。

**建議保留第 9 題作為研究框架，先將這輪成果視為後續比較的基礎，不把「R-QAOA 救回二分圖」當作已成立的創新。**

## 原本想驗證什麼

參考 Team 9 的圖結構研究，進一步測試：標準 QAOA 表現較弱的情況，能否透過遞迴縮圖改善？改善是否只出現在特定圖，還是能穩定延伸到其他結構？

R-QAOA 的方式是先估計兩個點偏向同組或異組，固定其中一個關係後縮小問題，反覆執行，到剩餘問題夠小時用古典方法解，再還原完整分組。

## 這輪已完成的工作

| 工作 | 範圍 |
|---|---|
| 八類圖比較 | 12、16 點，共 88 個案例 |
| 控制圖結構實驗 | 固定 16 點、32 邊、每點度數 4，20 組種子各 6 種改接程度，共 120 個案例 |
| 方法比較 | 標準 QAOA、三種 R-QAOA 設定、古典局部搜尋及 SDP／GW |
| 深度與雜訊補測 | 標準 QAOA 的 p=1、2、3；三張小圖的自訂 Aer 閘雜訊測試 |
| 結果驗證 | 主圖最佳值枚舉、原圖重新計分、解析公式與 Qiskit 核對 |

本輪主比較共 208 個案例；深度與雜訊是同圖上的補測，不能另算成大量獨立圖。全部使用本機模擬，沒有跑量子真機。原始圖由程式生成，並非 Team 9 原始圖的逐項重現。
'''),('結果與限制：目前沒有足夠的突破證據','''
## 1. R-QAOA 可以解對很多圖，但古典對照更強

| 方法 | 找到最佳解的案例數／208 |
|---|---:|
| R-QAOA | 197 |
| 古典局部搜尋，32 個起點 | 206 |
| SDP／Goemans–Williamson，256 次取整 | 207 |

這是指定樣本的結果，不是所有圖的勝率。三種 R-QAOA 設定在這批圖的最終分數相同，受限參數域版本沒有增加成功案例。各方法總資源尚未嚴格匹配，因此不把此表當公平成本下的最終排名；它已足以說明，目前不能主張 R-QAOA 有量子優勢。

## 2. 對平均分數的改善，不能直接當成找答案的優勢

標準 QAOA 量測 256 次、保留最好答案後，在這批小圖上本來就常接近最佳解。因此若只把「原 QAOA 平均分數」和「R-QAOA 最後那一個答案」相比，改善幅度容易被放大。

控制實驗還看到：改接連線後，QAOA 平均近似比由約 0.662 升至 0.804，但實際切開的邊比例由約 66.2% 降到 64.2%。比值變好的一個重要原因，是圖的最佳切割值也降低了。這是值得保留的評估發現，但目前不是新演算法的突破。

## 3. 失敗原因已找到，但改善尚未驗證

R-QAOA 可能在前幾步固定錯關係，直接排除最佳答案。例如一張最佳值為 25 的圖，第一步固定異組後，最高只剩 23，最後也只得到 23。這解釋了部分失敗，尚未產生能穩定修正它的新方法。

## 4. 幾個直覺改法已有前人研究

R-QAOA 解二分圖已有研究；用強化學習調整決策、自適應分配量測、加入回溯或分支搜尋，也都有相關工作。這不表示不能沿用，但不能只換一個組合或名稱，就當成我們的創新。

**目前真正缺少的是：一個明確的改動，在新的圖上、相同預算下，確實帶來可重現的收益，並能說清楚與既有研究的差別。**
'''),('下一步建議：先找兩條方法的結構互補性','''
## 保持第 9 題主線

團隊目前分別嘗試 R-QAOA 與 Warm-start QAOA。建議下一步先合併到同一套比較，回答：

> R-QAOA 和 Warm-start 各自適合哪些圖？能否找到一種有解釋力的結構條件，讓我們知道什麼時候應該用哪個方法？

這是下一個待驗證問題，不是已證明的新方法，也不是保證尚無人研究。即使後續找到互補性，「方法選擇」仍需再與相關文獻比對。

## 具體執行順序

**第一步：同圖比較。** 兩條線共用現有邊清單與最佳值，統一回報平均分數、固定量測數的最好答案、成功率與成本。不要讓不同圖、不同深度或不同抽樣數混在一起比較。

**第二步：看是否真的互補。** 找出「R-QAOA 成功、Warm-start 失敗」及反向情況，檢查它們是否對應可描述的圖結構，而不是少數偶然案例。先用簡單結構指標，不另開大型機器學習題目。

**第三步：只挑一個有證據的缺口改。** 如果某類結構反覆讓同一方法失敗，再針對該問題選擇改良；回溯、加深或初始化都可作候選，但須承認既有來源、納入成本，不能先當原創。

**第四步：鎖定規則後測新圖。** 現有案例只用於開發。新圖不能用來反覆調規則；同時比較原方法、簡單選擇規則、直接古典求解及同預算重啟，確認改善不是因為多試幾次。

## 什麼結果才值得繼續深化？

若能在未見圖上，依結構選擇或改良方法後，穩定優於固定使用單一量子方法，且成本可接受，就有值得深入的實驗結果。若要進一步宣稱量子優勢，還必須超過有競爭力的古典對照。

若兩條線沒有明顯互補，或改善只來自額外運算，就不繼續把 R-QAOA 包裝成突破。這輪資料仍可保留作第 9 題的基準與負結果，但是否足以主打參賽，需要重新評估。

## 本次希望團隊帶走的結論

**R-QAOA 路線已跑通；目前尚未形成足夠突破。下一步最有用的工作，是把兩條線放在相同條件下比較，先找出確實存在的結構缺口，再決定是否投入改良。**
'''),('補充：方法範圍與參考資料','''
## 閱讀數據時要保留的限制

主比較的 R-QAOA 固定 p=1，逐次消去一個變數，到 4 點精確求解。p=2、3 補測僅針對標準 QAOA；雜訊補測使用理想參數訓練後的自訂 Aer 執行噪音，不是 IBM 校準硬體，也沒有端到端含噪訓練。

受限參數域版本參考既有文獻，但本程式採數值搜尋，沒有全域參數最優證書。SDP／GW 對照也是數值求解與隨機取整，不把 SDP 目標值當嚴格最佳性證書。

控制圖雖固定點數、邊數與度數，改接時仍同時影響三角形、最佳值與其他結構，不能將結果全部歸因於單一指標。

完整技術報告、圖表、程式與原始 CSV 已整理，可供後續對照；本文件只保留團隊討論需要的結果與建議，尚未包含 Warm-start 的實測結果。

## 主要文獻

[Herrman 等：Impact of Graph Structures for QAOA on MaxCut](https://arxiv.org/abs/2102.05997)。提供第 9 題圖結構比較的研究背景。

[Bravyi 等：Obstacles to Variational Quantum Optimization](https://arxiv.org/abs/1910.08980)。R-QAOA 的原始研究。

[Bae 等：Modified Recursive QAOA for Exact Max-Cut Solutions on Bipartite Graphs](https://arxiv.org/abs/2408.13207)。二分圖與受限參數域方法。

[Patel 等：Reinforcement Learning Assisted Recursive QAOA](https://arxiv.org/abs/2207.06294)。失敗分析與學習式改良。

[Finžgar 等：Quantum-Informed Recursive Optimization Algorithms](https://arxiv.org/abs/2308.13607)。量子資訊引導縮圖與回溯。

[Lee、Kim：Adaptive Shot Allocation for Recursive QAOA](https://arxiv.org/abs/2605.26544)。逐步量測預算分配的相關預印本。

[Hadfield 等：Iterative Quantum Optimization: From Greedy Reductions to Branch-and-Bound](https://arxiv.org/abs/2611.00678)。可逆搜尋與界限方法的新預印本；其主張未於本次實驗重現。
''')]

def build():
    OUT.mkdir(exist_ok=True);PDF.mkdir(exist_ok=True)
    md=OUT/'第9題_RQAOA研究回報_團隊版.md'
    md.write_text('# 第 9 題｜Recursive QAOA 初步研究回報\n\n'+'\n\n---\n\n'.join('## '+t+'\n\n'+b.strip() for t,b in pages),encoding='utf-8')
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    from reportlab.platypus import SimpleDocTemplate,Paragraph,Table,TableStyle,PageBreak,Spacer
    from reportlab.lib.styles import ParagraphStyle
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4
    pdfmetrics.registerFont(TTFont('TeamRegular',r'C:\Windows\Fonts\msjh.ttc',subfontIndex=0));pdfmetrics.registerFont(TTFont('TeamBold',r'C:\Windows\Fonts\msjhbd.ttc',subfontIndex=0));pdfmetrics.registerFontFamily('TeamRegular',normal='TeamRegular',bold='TeamBold')
    body=ParagraphStyle('Body',fontName='TeamRegular',fontSize=10.8,leading=17.3,spaceAfter=10,wordWrap='CJK',textColor=colors.HexColor('#24334a'))
    title=ParagraphStyle('Title',parent=body,fontName='TeamBold',fontSize=19,leading=27,spaceAfter=20,textColor=colors.HexColor('#135f66'))
    heading=ParagraphStyle('Head',parent=body,fontName='TeamBold',fontSize=12.7,leading=18,spaceBefore=6,spaceAfter=8)
    small=ParagraphStyle('Small',parent=body,fontSize=9,leading=13,spaceAfter=3)
    def p(text,style=body):
        text=html.escape(text);text=re.sub(r'\*\*(.*?)\*\*',r'<b>\1</b>',text);text=re.sub(r'\[([^\]]+)\]\((https?://[^)]+)\)',r'<link href="\2" color="#167d82">\1</link>',text);return Paragraph(text,style)
    flow=[]
    for idx,(t,b) in enumerate(pages):
        if idx:flow.append(PageBreak())
        flow.append(p(t,title))
        for block in b.strip().split('\n\n'):
            if block.startswith('## '):flow.append(p(block[3:],heading))
            elif block.startswith('|'):
                lines=[l for l in block.splitlines() if not re.match(r'^\|[\s:|\-]+\|$',l)]
                data=[[p(v.strip(),small) for v in l.strip('|').split('|')] for l in lines]
                table=Table(data,colWidths=[155,335],hAlign='LEFT',repeatRows=1);table.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e6f0f0')),('VALIGN',(0,0),(-1,-1),'TOP'),('LINEBELOW',(0,0),(-1,0),.7,colors.HexColor('#167d82')),('LINEBELOW',(0,1),(-1,-1),.3,colors.HexColor('#dce4e9')),('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6),('LEFTPADDING',(0,0),(-1,-1),6),('RIGHTPADDING',(0,0),(-1,-1),6)]));flow.extend([table,Spacer(1,11)])
            else:flow.append(p(block.lstrip('> ').replace('\n',' ')))
    def footer(canvas,doc):
        canvas.setFont('TeamRegular',8);canvas.setFillColor(colors.HexColor('#6b788a'));canvas.drawString(48,29,'Johnny｜第 9 題研究回報｜2026-11-07');canvas.drawRightString(547,29,str(doc.page));canvas.setStrokeColor(colors.HexColor('#dce4e9'));canvas.line(48,43,547,43)
    pdf=PDF/'第9題_RQAOA研究回報_團隊版.pdf';doc=SimpleDocTemplate(str(pdf),pagesize=A4,leftMargin=48,rightMargin=48,topMargin=45,bottomMargin=58,title='第9題 Recursive QAOA 初步研究回報',author='Johnny');doc.build(flow,onFirstPage=footer,onLaterPages=footer)
    print(md);print(pdf)

if __name__=='__main__':build()
