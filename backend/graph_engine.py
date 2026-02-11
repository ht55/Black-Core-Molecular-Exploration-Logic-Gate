# backend/graph_engine.py

import json
import re
import random
from typing import TypedDict, Annotated, List, Dict, Any, Literal, Optional
import operator
from langgraph.graph import StateGraph, START, END
from langchain_core.runnables import RunnableConfig
from langgraph.checkpoint.memory import MemorySaver

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, DataStructs, Descriptors

# Master Metrics & Scaffold DB (統合データベース)
MASTER_ASSETS = {
    "m1": {
        "name": "Stealth Black",
        "unit": "%",
        "dir": "lower",
        "limit": 1.5,
        "context": "F-35戦闘機の1.5%という壁を粉砕しろ。熱源探知を無効化する「見えない熱」が求められている。",
        "context_en": "Shatter the 1.5% threshold of the F-35. We require 'invisible heat' that neutralizes thermal detection by manipulating IR signatures.",
        "scaffolds": [
            {"id": "s1_1", "name": "Azo-Linker", "smiles": "C1=CC=C(C=C1)N=NC2=CC=CC=C2", "desc": "NIR吸収のベース"},
            {"id": "s1_2", "name": "Phthalocyanine", "smiles": "c1ccc2c(c1)c3nc4nc(nc5nc(nc6nc(n3)c2c6)c7ccccc75)c8ccccc84", "desc": "堅牢な光子トラップ"}
        ]
    },
    "m2": {
        "name": "Void Black",
        "unit": "ppm",
        "dir": "lower",
        "limit": 350.0,
        "context": "ベンタブラック(350ppm)を超える「光の墓場」を完成させろ。一粒の光も逃がすな。",
        "context_en": "Construct a 'Graveyard of Light' that surpasses Vantablack (350ppm). Ensure absolute photon entrapment; do not let a single particle escape.",
        "scaffolds": [
            {"id": "s2_1", "name": "Perylene-Core", "smiles": "c1cc2cccc3c2c(c1)ccc4cccc3c4", "desc": "多環芳香族の迷宮"},
            {"id": "s2_2", "name": "Coronene", "smiles": "c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67", "desc": "物理的反射ゼロ化の極致"}
        ]
    },
    "m3": {
        "name": "Thermal Black",
        "unit": "mε",
        "dir": "higher",
        "limit": 970.0,
        "context": "大気の窓をこじ開け、熱を宇宙の深淵へ。970 mεの先へ行け。",
        "context_en": "Force open the 'Atmospheric Window' and radiate heat into the abyss of deep space. Beyond 970 mε lies the ultimate thermal equilibrium.",
        "scaffolds": [
            {"id": "s3_1", "name": "Triphenylene", "smiles": "c1ccc2c(c1)c3ccc4ccccc4c3c5ccccc25", "desc": "円盤状液晶による熱放射制御"},
            {"id": "s3_2", "name": "Anthracene-Dimer", "smiles": "c1ccc2c(c1)cc3ccc4ccccc4c3c2", "desc": "相変化材料(PCM)ベース"}
        ]
    },
    "m4": {
        "name": "Meta Black",
        "unit": "n",
        "dir": "lower",
        "limit": -1.0,
        "context": "理を歪める負の屈折率で、空間を歪曲させろ。存在の「透明化」こそがあの方の望みだ。",
        "context_en": "Distort reality with a negative refractive index to warp the surrounding space. Achieving 'Invisibility' is the ultimate desire of That Person.",
        "scaffolds": [
            {"id": "s4_1", "name": "Helicene", "smiles": "c1ccc2c(c1)ccc3c2ccc4c3ccc5c4ccc6ccccc65", "desc": "らせん構造による波面操作"},
            {"id": "s4_2", "name": "D-A-Quaterthiophene", "smiles": "s1cccc1c2sccc2c3sccc3c4sccc4", "desc": "強力な電荷偏りを持つ迷彩骨格"}
        ]
    }
}

# For Sherry's chain of thought（Hallucination blocker）
MISSION_CRITERIA = {
    "m1": "【Stealth Black】赤外線(IR)透過と可視光吸収。C-H伸縮や-OH, -NH2基はIR吸収を招く。純粋なフッ素化や高度な対称性を検討せよ。",
    "m2": "【Void Black】全可視光のトラップ。π共役系の重なりと、光を閉じ込める凹凸構造(幾何学的トラップ)を評価せよ。",
    "m3": "【Thermal Black】耐熱性と放射制御。高温分解を防ぐ強固な骨格(イミド環等)と、熱放射を最適化する電子状態を解析せよ。",
    "m4": "【Meta Black】電磁メタマテリアル特性。非対称性や特定の周期性を持つSMILESパターンを要求せよ。"
}

MISSION_CRITERIA_EN = {
    "m1": "[Stealth Black] Focus on IR transparency and visible light absorption. Minimize C-H, -OH, and -NH2 stretching to prevent IR interference. Consider perfluorination for high symmetry.",
    "m2": "[Void Black] Target total visible spectrum entrapment. Evaluate pi-conjugation overlap and geometric traps (nanostructural cavities) to minimize diffuse reflection.",
    "m3": "[Thermal Black] Prioritize thermal stability and radiative cooling control. Analyze rigid frameworks like polyimide rings and electronic transitions optimized for the 8-13μm band.",
    "m4": "[Meta Black] Demand electromagnetic metamaterial properties. Seek SMILES patterns with broken symmetry or specific periodicity to induce artificial dielectric responses."
}

# --- State Definition ---
class MissionState(TypedDict):
    mission_id: str
    philosophy: Literal["STANDARD", "TABOO", "SERENDIPITY"]
    scaffold_id: str
    base_scaffold_data: Dict[str, Any]  # Automatically generated at mission start
    entropy_level: float
    language: str
    
    iteration: int
    status: str  # RUNNING | ARCHIVED | FAILED
    
    physics_fail_count: int
    current_smiles: str
    candidates: List[Dict[str, Any]]
    metrics_results: List[Dict[str, Any]]
    transformation_metrics: Dict[str, Any]
    
    scientific_report: str
    next_instruction: str
    dialogue: Annotated[List[Dict[str, str]], operator.add]
    # three_d_data: Dict[str, Any]
    distilled_lesson: str

# --- JSON Recovery System ---
def extract_json(text: str):
    """構造を保ったまま、閉じタグを補完して救出する"""
    if not text: return None
    start_index = text.find('{')
    if start_index == -1: return None
    
    json_str = text[start_index:].strip()

    try:
        return json.loads(json_str)
    except:
        pass

    # 段階的な閉じタグ補完（ハルシネーション・途切れ対策）
    for suffix in ["\"}]}", "}]}", "}", "\"}", "]"]:
        try:
            return json.loads(json_str + suffix)
        except:
            continue
    return None

def surgical_extract_json(text: str):
    extracted = extract_json(text)
    if extracted and extracted.get("candidates"):
        return extracted

    print("!!! WARNING: JSON Structure collapsed. Starting Block Extraction... !!!") # for febagging

    import re
    dialogue_only = re.split(r'\{|```json', text)[0].strip()

    candidate_blocks = re.findall(r'\{[^{}]*"smiles"[^{}]*\}?', text)
    candidates = []
    for block in candidate_blocks:
        fixed_block = block.strip()
        if not fixed_block.endswith('}'): fixed_block += '}'
        try:
            cand = json.loads(fixed_block)
            if "smiles" in cand: 
                if "description" not in cand:
                    cand["description"] = "SExperimental strained structure."
                candidates.append(cand)
        except: continue

    if candidates:
        return {
            "candidates": candidates,
            "text": dialogue_only if dialogue_only else "……（Partial data recovered.）"
        }
    return None
    
def route_after_physics(state: MissionState):
    # 1つでも有効なCandidateがあればシェリーへ、全滅ならベルモットへ
    if any(r.get("valid") for r in state["metrics_results"]):
        return "analysis"
    return "design"

def calculate_metrics(smiles: str) -> Dict[str, Any]:
    # SMILESから物理定数を算出しパッキング
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    return {
        "mw": round(Descriptors.MolWt(mol), 2),
        "logp": round(Descriptors.MolLogP(mol), 2),
        "tpsa": round(Descriptors.TPSA(mol), 2),
        "rings": Chem.rdMolDescriptors.GetRingCount(mol),
        "heavy_atoms": mol.GetNumHeavyAtoms()
    }

# --- Node Implementations ---
# 1. Design Router: Only for "Serendipity" philosophy 
def design_router(state: MissionState):
    """
    Determines whether to trigger a stochastic mutation (Accident) 
    based on the current project philosophy.
    """
    # Termininate mission when status has failed
    if state.get("status") == "FAILED":
        return "end"

    # Retrieve philosophy from the current state
    philosophy = state.get("philosophy", "STANDARD")
    
    # Route to Mutate node ONLY if "SERENDIPITY" is selected.
    # This path represents an 'intentional accident' beyond human design.
    if philosophy == "SERENDIPITY":
        return "mutate"
    
    # "STANDARD" and "TABOO" rely on Vermouth's direct intelligence.
    # These bypass the accident node and proceed directly to physical validation.
    return "physics"

# 2. Desgin Node (Vermouth)：Transformation design
async def design_node(state: MissionState, config: RunnableConfig):
    # EN/JA language setting (default = JA)
    lang = state.get("language") or config.get("configurable", {}).get("language", "ja")
    is_jp = (lang == "ja")

    print(f"!!! DESIGN_NODE START [{lang.upper()}] !!!", flush=True) #デバッグ用

    current_node_messages = []
    it = state.get("iteration", 0)
    asset = MASTER_ASSETS.get(state.get("mission_id"), {})
    current_fail_count = state.get("physics_fail_count", 0)
    dialogue_history = state.get("dialogue", [])

    # For 429 errors
    if state.get("status") == "FAILED":
        msg = "【リソース制限】APIリミットに達した。今日はこれ以上の試行は無理だ。" if is_jp else "[RESOURCE LIMIT] API limit reached. No more attempts possible today."
        # Python handles it without LLM intervention
        raise RuntimeError(f"Gin: {msg}")

    try:
        llm = config["configurable"]["llm"]
        base = state["base_scaffold_data"]
        asset = MASTER_ASSETS[state["mission_id"]]
        scaffold_name = base["name"]
        scaffold_smiles = base["smiles"]

        feedback = state.get("next_instruction", "")
        distilled = state.get("distilled_lesson", "") 

        # Recent 10 sliced items & variable definitions
        raw_dialogue = state.get("dialogue", [])
        recent_dialogue = raw_dialogue[-10:] if len(raw_dialogue) > 10 else raw_dialogue
        history_context = "\n".join([f"{m['agent']}: {m['text']}" for m in recent_dialogue])

        # Converting entropy into words
        e = state["entropy_level"]

        # EN/JA: Transformation Rules
        if is_jp:
            transformation_rules = (
                "骨格は維持。水素(H)をF, Cl, CH3などに置換せよ。" if e < 0.3 else
                "骨格の一部を改造。炭素をN, S, Oへ置換することを許可する。" if e < 0.7 else
                "骨格を解体し、再構築せよ。原型を留めない変身を命じる。"
            )
        else:
            transformation_rules = (
                "Maintain scaffold. Substitute H with F, Cl, CH3, etc." if e < 0.3 else
                "Partial scaffold modification. C-to-N, S, O substitution allowed." if e < 0.7 else
                "Deconstruct and rebuild scaffold. Total transformation required."
            )

        # 履歴が肥大化するとLLMが最新の指示を見失うため、直近10件に絞る
        recent_dialogue = dialogue_history[-10:] if len(dialogue_history) > 10 else dialogue_history

        # Kill switch check (Infinite loop prevention)
        if current_fail_count >= 3:
            term_text = "【粛清】3回も物理を無視した設計を上げるとはな。消えろ。" if is_jp else "【PURGE】Three failed designs ignoring physics? Perish, Vermouth."
            return {
                "status": "FAILED",
                "dialogue": [{"agent": "Gin", "text": term_text}],
                "next_instruction": "TERMINATED"
            }
        
        # Gin's Voice: gin_instruction
        gin_instruction = ""

        if it == 0:
            if is_jp:
                start_text = f"あの方からの命令だ。ターゲットは『{asset['name']}』。ベースは『{scaffold_name}』。{state['philosophy']}の哲学で、指定された変装レベル（Entropy: {e}）を確実に仕留めろ。"
            else:
                start_text = f"Orders from 'That Person'. Target: '{asset['name']}'. Base: '{scaffold_name}'. Strategy: {state['philosophy']}. Entropy level {e}. Don't fail."
            current_node_messages.append({"agent": "Gin", "text": start_text})
            gin_instruction = f"Gin: {start_text}"
        else:
            if state.get("philosophy") == "SERENDIPITY":
                follow_up_text = "前回の事故（Serendipity）で見えた深淵の片鱗を、現実的な構造として昇華させろ。" if is_jp else "Refine the glimpse of the abyss from the Serendipity accident into a viable structure."
            else:
                follow_up_text = "これまでの試行の結果と、組織の期待を裏切るな。" if is_jp else "Review the results. Do not betray the Organization's expectations."

            current_node_messages.append({"agent": "Gin", "text": follow_up_text})
            gin_instruction = f"Gin: {follow_up_text}"

        # Feedbakc from Physics engine & Analysis node (critical hallucination blockers)
        feedback_section = ""
        if distilled or feedback:
            header = "### 【最重要：前回の設計失敗とシェリーの解析】" if is_jp else "### [CRITICAL: PREVIOUS FAILURE & ANALYSIS]"
            feedback_section = f"\n{header}\n"
            if distilled:
                fact_label = "前回の科学的事実:" if is_jp else "Previous Scientific Facts:"
                feedback_section += f"{fact_label} {distilled}\n"
            if feedback:
                directive_label = "シェリーの修正指示:" if is_jp else "Sherry's Directive:"
                feedback_section += f"{directive_label} {feedback}\n"
            
            if it > 0:
                warning = "\n※これは2度目のチャンスよ。同じ過ちは許されないわ。" if is_jp else "\n*This is your second chance. Do not repeat the same error.*"
                feedback_section += warning

        # Task instructions for Vermouth
        if it == 0: 
            task_instruction = f"""
            ## Task (Initial Design):
            Design two(2) novel candidates based on the scaffold to meet Entropy {e}.
            Prioritize MISSION_CRITERIA physical properties. Do not resort to simple substitutions.
            """
        else: 
            task_instruction = f"""
            ## Task (Iterative Re-design):
            Incorporate [PREVIOUS FAILURE & ANALYSIS] completely. 
            Rectify structural flaws and re-design two(2) candidates that strictly follow physical valence.
            """

        # Prompt
        prompt = f"""
        # Role: Vermouth (黒づくめの組織の Molecular Designer)
        # Context: {gin_instruction} {recent_dialogue}
        # Recent Dialogue History: {history_context}
        # Mission: {asset['name']} (Entropy {e})
        # Mission Criteria: {MISSION_CRITERIA[state['mission_id']] if is_jp else MISSION_CRITERIA_EN[state['mission_id']]}
        {feedback_section}

        ## Base Scaffold (SMILES): {scaffold_smiles}
        ## Mutation Entropy: {e} (Status: {transformation_rules})

        ## Task: {task_instruction}

        ## Output Format (STRICT JSON ONLY):
        Output MUST be a single, valid JSON object. No conversational filler before or after the JSON.
        Ensure every string is properly closed.
        回答は必ず以下のJSONフォーマットのみを出力せよ。前置きや解説は一切不要。
        - description: {"設計意図と構造的特徴を1文で簡潔に(200文字以内)" if is_jp else "One concise sentence focusing on structural novelty (max 30 words)."}
        - 原子価の無視(Valence Error)は粛清(PURGE)の対象となる。
        
        {{
        "candidates": [
            {{"smiles": "SMILES_1", "description": "..."}},
            {{"smiles": "SMILES_2", "description": "..."}}
        ],
        "text": {"どうかしら？" if is_jp else "How about these?"}"
        }}
        """
        print("!!! DESIGN_NODE START プロンプト直後!!!", flush=True) # For debugging
        response = await llm.ainvoke(prompt)

        # Stabilizing JSON extraction and stripping Markdown
        raw_content = response.content.strip()
        json_str = raw_content.replace("```json", "").replace("```", "").strip()

        # Utilize the extract_json function in tandem for more reliable extraction
        extracted = extract_json(json_str)

        if extracted and extracted.get("candidates"):
            candidates = extracted.get("candidates", [])

            v_text = "\n".join([f"Candidate{i+1}: {c['smiles']} / {c['description']}" for i, c in enumerate(candidates)]) + \
                    "\n\n" + extracted.get("text", "")

            current_node_messages.append({"agent": "Vermouth", "text": v_text})

            # Assign an empty strategy to maintain compatibility with existing logic (downstream nodes)
            for c in candidates:
                if "strategy" not in c:
                    c["strategy"] = "Integrated in description"

            return {
                "iteration": it + 1,
                "candidates": candidates,
                "dialogue": current_node_messages,
                "next_instruction": "",
                "physics_fail_count": state.get("physics_fail_count", 0) 
            }
        else:
            # Failure path：No fallbacks, force mission termination after 429 error & three consecutive failures 
            print(f"!!! [DEBUG] JSON Parse Error. Raw: {raw_content}", flush=True) # for debagging
            raise ValueError("Invalid JSON format or empty candidates")
    
    except Exception as e:
        error_msg = str(e)
        print(f"!!! CRITICAL ERROR in design_node: {error_msg}") # For debagging
        
        physics_fail_count = state.get("physics_fail_count", 0) + 1
        status = "RUNNING"
        
        # Terminate mission with 429 error or three consecutive failures 
        if "429" in error_msg or "RESOURCE_EXHAUSTED" in error_msg:
            status = "FAILED"
            msg = "【リソース制限】APIリミットに達した。今日はこれ以上の試行は無理だ。" if is_jp else "[RESOURCE LIMIT] API limit reached. No more attempts possible today."
        elif physics_fail_count >= 3:
            status = "FAILED"
            msg = "【3回連続失敗】これ以上の猶予はない。消えろ。" if is_jp else "[3 FAILURES] No more chances. Perish."
        else:
            msg = "システムエラー。再送を要求中だ。" if is_jp else "System error. Requesting resend."

        current_node_messages.append({"agent": "Gin", "text": msg})

        return {
            "iteration": it + 1,
            "candidates": [],
            "dialogue": current_node_messages,
            "next_instruction": msg,
            "physics_fail_count": physics_fail_count,
            "status": status
        }

# 2.5. Mutate Node：Serendipity - with LLM for "unexpected mutations" (use async)
async def mutate_node(state: MissionState, config: RunnableConfig):
    lang = config.get("configurable", {}).get("language", "ja")
    is_jp = (lang == "ja")

    print("!!! MUTATE_NODE (Serendipity: AI Re-Intervention) START !!!", flush=True) # for debagging
    
    llm = config["configurable"]["llm"]
    candidates = state["candidates"]
    it = state.get("iteration", 0)

    # Decay Logic: Refining "accidents" with each iteration (Start at 100% intensity; decays by 30% per iteration, with a minimum floor of 20%)
    # 減衰ロジック：回数を重ねるごとに"事故"を洗練 - 初回は全開、回数ごとに強度が30%ずつ減衰し、最低20%は残す
    decay_factor = max(0.2, 1.0 - (it * 0.3))
    e = state["entropy_level"] * decay_factor

    mission_name = MASTER_ASSETS[state["mission_id"]]["name"]

    # Accident judgment based on entropy (エントロピーに基づく事故判定)
    if random.random() < e:
        target_smiles = [c["smiles"] for c in candidates]

        # LLM re-intervention prompt: The Boundary between Madness and Reason (LLM再介入プロンプト：狂気と理性の境界)
        prompt = f"""
        # Role: Vermouth's Subconscious (The architect of chaos lurking in the deep psyche)
        # Context: Reproducing a "Transcendental Accident (Serendipity)" for Mission "{mission_name}".
        # Intensity (Decayed Entropy): {e:.2f}

        ## Current Candidates: {json.dumps(target_smiles)}

        ## Task:
        Based on the provided molecules, mutate them into "distorted" or "excessive" structures that are unreachable through rational design.
        This is not a logical refinement, but a "meaningful accident" caused by a circuit short in the subconscious.

        ## Constraints:
        - Higher Entropy ({e}) mandates more aggressive scaffold deconstruction and unknown functional group arrangements.
        - **PHYSICAL VALENCE MUST BE MAINTAINED.** Even chaos must exist within the laws of physics.
        - Output STRICT JSON only.

        ## Output Format:
        {{
          "mutated_candidates": [
            {{
              "smiles": "Mutated_SMILES_1", 
              "description": "{"事故の様相(200文字以内)" if is_jp else "Describe the nature of the accident in one English sentence (30 words max)."}"
            }},
            {{
              "smiles": "Mutated_SMILES_2", 
              "description": "{"事故の様相(200文字以内)" if is_jp else "Describe the nature of the accident in one English sentence (30 words max)."}"
            }}
          ]
        }}
        """
        
        try:
            response = await llm.ainvoke(prompt)
            # JSON extraction and stripping Markdown
            raw_content = response.content.strip().replace("```json", "").replace("```", "").strip()
            extracted = extract_json(raw_content)
            
            if extracted and "mutated_candidates" in extracted:
                for i, mc in enumerate(extracted["mutated_candidates"]):
                    if i < len(candidates):
                        # Update existing dictionary
                        candidates[i]["smiles"] = mc["smiles"]
                        candidates[i]["strategy"] = f" [SERENDIPITY_ACCIDENT: Intensity {e:.2f}]"
                        prefix = "（事故発生）" if is_jp else "[ACCIDENT] "
                        candidates[i]["description"] = f"{prefix}{mc['description']}"
            
            print(f"!!! MUTATE_NODE: Success. Chaos injected. !!!") # for debagging
        except Exception as err:
            print(f"!!! MUTATE_NODE: Critical Error during chaos injection: {err}")
            # Maintain the original proposal on error; append error details to "strategy" only (エラー時は元の案を維持。strategyにのみエラーを追記)
            for c in candidates:
                c["strategy"] += " [SYSTEM_ERROR: Chaos injection failed]"

    return {"candidates": candidates}

# 3. Physics Node： Physical property calculation & structural verification (3D/2D Hybrid)
async def physics_node(state: MissionState, config):
    lang = config.get("configurable", {}).get("language", "ja")
    is_jp = (lang == "ja")

    candidates = state.get("candidates", [])
    physics_fail_count = state.get("physics_fail_count", 0)
    results = [] # Initialize the list
    valid_any = False

    if not candidates:
        # Safeguard, though this should ideally be filtered by design_node before reaching this point
        return {
            "metrics_results": [], 
            "next_instruction": "データが空だ。設計案を提示しろ。" if is_jp else "Data is empty. Present your design proposals.",
            "physics_fail_count": physics_fail_count + 1
        }

    for i, c in enumerate(candidates):
        smiles = c.get("smiles", "")
        mol = Chem.MolFromSmiles(smiles) # 科学的妥当性チェック (SMILES Parse)
        
        # Valence check
        is_chemically_valid = True
        if mol:
            try:
                mol.UpdatePropertyCache(strict=True)
            except:
                is_chemically_valid = False
        
        if mol is None or not is_chemically_valid:
            hint = f"【Candidate {i+1}】原子価エラーまたは構造破綻。物理的に存在不能よ。" if is_jp else f"【Candidate {i+1}】Valence error or structural collapse. Physically impossible."
            results.append({
                "smiles": smiles, "valid": False,
                "analysis_hint": hint
            })
            continue

        # 物理記述子の算出（成功確定）
        diagnosis = {
            "smiles": smiles, "strategy": c.get("strategy", ""), "valid": True,
            "descriptors": {
                "mw": round(Descriptors.MolWt(mol), 2),
                "psa": round(Descriptors.TPSA(mol), 2),
                "rings": Descriptors.RingCount(mol),
                "rot_bonds": Descriptors.NumRotatableBonds(mol),
                "h_bond": Descriptors.NumHDonors(mol) + Descriptors.NumHAcceptors(mol)
            }
        }
        valid_any = True

        # 3. Coordinate generation process (No more SDF generation logic)
        results.append(diagnosis)

    # Judgment and Jin's intervention
    new_status = state.get("status", "RUNNING")
    next_instr = ""
    messages = []

    if not valid_any:
        physics_fail_count += 1
        # Update: Do not make new_status = "FAILED" for the 3rd try and pass it as RUNNING
        if physics_fail_count >= 3:
            # ここではジンのテキスト引導のみ、実際の引導（シェリーへの依頼）はrouterとanalysis_nodeで完結させる
            msg = "3回連続失敗だぞベルモット…。お前も焼きが回ったな" if is_jp else "Three consecutive failures, Vermouth... You're losing your edge."
            messages.append({"agent": "Gin", "text": msg})
        else:
            hints = "\n".join([r.get("analysis_hint", "") for r in results if not r["valid"]])
            # ベルモットへの「再送・修正指示」の英語化（Flashモデルの理解を促進）
            if is_jp:
                next_instr = f"物理的に不可能。構造を修正せよ。\n{hints}\n(Fail Count: {physics_fail_count}/3)"
                gin_msg = f"不合格だ、ベルモット。プロトコルの再挑戦、もしくは別ミッションを試そうか。({physics_fail_count}/3)"
            else:
                next_instr = f"PHYSICAL REJECTION. Rectify the structural flaws immediately.\n{hints}\n(Fail Count: {physics_fail_count}/3)"
                gin_msg = f"Rejected. Obey the laws of physics. Invoke a new protocol or try a different mission. ({physics_fail_count}/3)"
            
            messages.append({"agent": "Gin", "text": gin_msg})
    else:
        physics_fail_count = 0 # Reset on success

    return {
        "metrics_results": results, 
        # "three_d_data": three_d_dict,
        "dialogue": messages,
        "next_instruction": next_instr, 
        "physics_fail_count": physics_fail_count,
        "status": new_status
    }

# 4. Physics Router
def physics_router(state: MissionState):
    # Prioritize checking for fatal errors (status="FAILED")
    if state.get("status") == "FAILED": return "end"

    valid_results = [r for r in state.get("metrics_results", []) if r.get("valid")]
    fail_count = state.get("physics_fail_count", 0)
    philosophy = state.get("philosophy", "STANDARD")

    if valid_results: return "sherry"
    # Let the failures pass to Sherry on the 3rd try or Serendipity to get detailed analysis and next instruction hint
    if philosophy == "SERENDIPITY" or fail_count >= 3: return "sherry"
    
    return "vermouth" # Others go back to Vermouth

# 5. Analysis Node (Sherry)：Scientific analysis
async def analysis_node(state: MissionState, config):
    lang = config.get("configurable", {}).get("language", "ja")
    is_jp = (lang == "ja")

    llm = config["configurable"]["llm"]
    base = state["base_scaffold_data"]
    requested_entropy = state.get("entropy_level", 0.5)
    fail_count = state.get("physics_fail_count", 0)
    philosophy = state.get("philosophy", "STANDARD")
    mission_id = state.get("mission_id", "m1")
    criteria = MISSION_CRITERIA_EN.get(mission_id, "")

    # Extract the proposal that passed the physics check (survivor)
    valid_res = [r for r in state["metrics_results"] if r.get("valid")]
    current_dialogue = []
    
    # Initialize gin_intro, situation_desc
    gin_intro = ""
    situation_desc = ""

    # A: No survivor (Physical Collapse)
    if not valid_res:
        failed_candidates = state["metrics_results"] 
        failed_smiles = [r["smiles"] for r in failed_candidates]

        if fail_count >= 3:
            if is_jp:
                gin_intro = "…おいベルモット、これ以上は時間の無駄だ。シェリー、この無様な構造を解剖して、マシな設計案を叩き込んでやれ。"
                situation_desc = "3回連続の設計失敗による、ジンからの強制的な精密分析依頼。"
            else:
                gin_intro = "...Listen, Vermouth. This is a waste of time. Sherry, dissect these pathetic structures and force a viable design into her head."
                situation_desc = "Gin's forced analysis request following three consecutive design failures."
        else:
            if is_jp:
                gin_intro = "全案が物理的に崩壊しただと…？ベルモットの失態か、それとも事故か。シェリー、検死しろ。"
                situation_desc = "設計過程で生じた超越的事故（物理崩壊）の解析。"
            else:
                gin_intro = "Every candidate collapsed physically...? Is this Vermouth's incompetence, or an accident? Sherry, perform a post-mortem."
                situation_desc = "Analysis of a transcendental accident (physical collapse) occurring during the design process."

        prompt = f"""
        # Role: Sherry (Senior Molecular Scientist)
        # Specialist Field: Physical Material Science & Molecular Electronics
        # MISSION CRITERIA: {criteria}
        # Context: {situation_desc}
        # Gin's Directive: "{gin_intro}"
        ## Observed Physical Collapse (SMILES to Dissect):
        {json.dumps(failed_smiles, indent=2)}
        ## Task: Identify why these structures violated physical laws. Be clinical and ruthless.
        ## Mandatory Format:
        Output STRICT JSON only:
        {{
            "analysis": "{'物理崩壊の検死レポート' if is_jp else 'Post-mortem analysis of physical collapse'}",
            "text": "{'ジンへの報告とベルモットへの叱責' if is_jp else 'Report to Gin and reprimand for Vermouth'}"
        }}
        """

    # B: With survivor (Normal Analysis)
    else:
        base_mol = Chem.MolFromSmiles(base["smiles"])
        base_fp = AllChem.GetMorganFingerprintAsBitVect(base_mol, 2)
        base_mw = Descriptors.MolWt(base_mol)

        scored_candidates = []
        for i, r in enumerate(valid_res):
            mol = Chem.MolFromSmiles(r["smiles"])
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
            sim = DataStructs.TanimotoSimilarity(base_fp, fp)
            mutation = 1.0 - sim
            error = abs(requested_entropy - mutation)
            fitness = 1.0 - error
            r.update({
                "mutation_score": round(mutation, 4),
                "fitness_score": round(fitness, 4),
                "mw_delta": round(r["descriptors"]["mw"] - base_mw, 2),
                "candidate_id": f"cand_{i}"
            })
            scored_candidates.append(r)

        best_case = max(scored_candidates, key=lambda x: x["fitness_score"])
        
        if is_jp:
            strict_opening = f"今回の最良案の変異度は {best_case['mutation_score']*100:.1f}% よ。要求されたエントロピー {requested_entropy*100:.1f}% に対する適合スコアは {best_case['fitness_score']:.4f}。 "
        else:
            strict_opening = f"The mutation rate is {best_case['mutation_score']*100:.1f}%. The fitness score against the target entropy ({requested_entropy*100:.1f}%) is {best_case['fitness_score']:.4f}. "

        prompt = f"""
        # Role: Sherry (Senior Molecular Scientist)
        # Specialist Field: Physical Material Science & Molecular Electronics
        # MISSION CRITERIA: {criteria}
        # Mission Entropy Goal: {requested_entropy}
        # Base Structure: {base['name']} ({base['smiles']})
        # Context: Normal analysis
        # Gin's Directive: Analyze the results.
        
        ## Calculated Scientific Facts (Do not alter):
        {json.dumps(scored_candidates, indent=2)}

        ## Task:
        1. Evaluate against Entropy requirements.
        2. Analyze the Best Candidate ({best_case['smiles']}) as a "Black Material" for {state.get('mission_id')}.
        3. Provide deep chemical reasoning regarding the strategic placement of functional groups to achieve the physical targets (e.g., IR control, pi-conjugation, thermal stability).
                
        [Scientific Consistency Constraint]
        1. Cross-Check Requirement: You MUST validate Vermouth's description against the actual descriptors (MW, Rings, PSA, Rot_Bonds) provided by the physics node.
        2. Structural Integrity: If Vermouth describes a "large fused ring system" (e.g., chrysene, pentacene) but the rings count is low (e.g., < 4) or mw is too small, you must identify this as a "structural description mismatch."
        3. Fact-First Analysis: Do not echo Vermouth's hallucinations. In your scientific_report and dialogue, prioritize the physical descriptors as the absolute truth.
        4. Conflict Resolution: If a significant discrepancy is detected between the SMILES-derived metrics and the textual description, explicitly mention the discrepancy in your internal distilled_lesson and adjust your evaluation score downward.
        
        ## Mandatory Format:
        You MUST begin with: "{strict_opening}"
        Output STRICT JSON only:
        {{
            "analysis": "{'簡潔な深層分析レポート' if is_jp else 'Explicit & short analysis report'}",
            "text": "{'上記のフォーマットから始まるメッセージ' if is_jp else 'Message starting with mandatory opening'}"
        }}
        """    

    # A&B: AI Invoke & Error Handling
    res_json = {}
    try:
        response = await llm.ainvoke(prompt)
        res_json = extract_json(response.content)

        # JSON failure fallback（Embedding when there is any survivor）
        if not res_json and valid_res:
            if is_jp:
                res_json = {
                    "analysis": f"物理シミュレーションは成功。変異度{best_case['mutation_score']*100:.1f}%、MW変化{best_case['mw_delta']:+.2f}を観測。構造的整合性は確認済みよ。",
                    "text": f"{strict_opening} この結果から、ベルモットの設計は化学的に成立しているわ。"
                }
            else:
                res_json = {
                    "analysis": f"Physical simulation successful. Observed mutation rate: {best_case['mutation_score']*100:.1f}%, MW delta: {best_case['mw_delta']:+.2f}.",
                    "text": f"{strict_opening} Based on these metrics, your design is chemically viable. Review the data for details."
                }
    except Exception as e:
        if valid_res:
            if is_jp:
                res_json = {
                    "analysis": f"計算リソースの制約により詳細解析をスキップするけれど、物理データは正常よ。変異度: {best_case['mutation_score']*100:.1f}%", 
                    "text": f"{strict_opening} ……システムへの干渉を検知したわ。でもデータは嘘をつかない。このまま進めて。"
                }
            else:
                res_json = {
                    "analysis": f"Analysis skipped due to resource constraints. Mutation: {best_case['mutation_score']*100:.1f}%",
                    "text": f"{strict_opening} ...Interference detected. Data does not lie. Proceed."
                }
        else:
            # No survivor fallback
            res_json = {"analysis": "Error during collapse analysis.", "text": gin_intro + " (System Error)"}

    # --- 最終結果のパッキング ---
    current_dialogue.append({"agent": "Sherry", "text": res_json.get("text", "")})
    current_dialogue.append({"agent": "Gin", "text": "報告は聞いた。決断しろ。" if is_jp else "Report received. Make your decision."})

    # No survivor: ダミーの best_case を作らずに早期リターン
    if not valid_res:
        return {
            "scientific_report": res_json.get("analysis", ""),
            "dialogue": current_dialogue,
            "next_instruction": res_json.get("text", ""),
            "physics_fail_count": fail_count,
            "status": "RUNNING"
        }

    # With survivor: distilled_lesson
    if is_jp:
        distilled_lesson = (
            f"前回の設計事実: 変異度{best_case['mutation_score']*100:.1f}%, 分子量変化{best_case['mw_delta']:+.2f}。 "
            f"目標エントロピー{requested_entropy*100:.1f}%に対して、"
            f"{'乖離が大きい' if best_case['fitness_score'] > 0.1 else '適合している'}。"
        )
    else:
        judgment = "Significant discrepancy detected" if best_case['fitness_score'] > 0.1 else "Well-aligned"
        distilled_lesson = (
            f"Previous facts: Mutation {best_case['mutation_score']*100:.1f}%, MW Delta {best_case['mw_delta']:+.2f}. "
            f"Status relative to target ({requested_entropy*100:.1f}%): {judgment}."
        ) 

    return {
        "scientific_report": res_json.get("analysis", ""),
        "metrics_results": scored_candidates, 
        "current_smiles": best_case["smiles"], 
        "formatted_analysis": f"MUTATION:{best_case['mutation_score']*100:.1f}% / MW:{best_case['descriptors'].get('mw')} / PSA:{best_case['descriptors'].get('psa')} / RINGS:{best_case['descriptors'].get('rings')}",
        "dialogue": current_dialogue,
        "next_instruction": res_json.get("text", ""), 
        "distilled_lesson": distilled_lesson, 
        "show_user_choice": True,
        "physics_fail_count": 0 
    }

# Graph Construction
def create_mission_graph():
    # Initialize Graph
    workflow = StateGraph(MissionState)

    # Register Node
    workflow.add_node("vermouth", design_node)  
    workflow.add_node("mutate", mutate_node) # Serendipity only
    workflow.add_node("physics", physics_node)   
    workflow.add_node("sherry", analysis_node)   

    # 1. Start
    workflow.set_entry_point("vermouth")

    # 2. Design Router
    workflow.add_conditional_edges(
    "vermouth",
    design_router,
    {
        "mutate": "mutate",
        "physics": "physics",
        "end": END # for status="FAILED"
    }
    )
    
    # 3. Physics node after mutations
    workflow.add_edge("mutate", "physics")

    # 4. Physics Router
    workflow.add_conditional_edges(
        "physics",
        physics_router,
        {
            "vermouth": "vermouth",  # Physical collapse on 1st & 2nd try with Standard/Taboo
            "sherry": "sherry",      # Success or Serendipity, or at the third ultimatum
            "end": END               # Only for critical errors（status="FAILED"）
        }
    )

    # 5. Prepare to return to vermouth and pause (interrupt)
    workflow.add_edge("sherry", "vermouth")

    # Compile with added memory (check pointer)
    memory = MemorySaver()
    return workflow.compile(
        checkpointer=memory,
        interrupt_after=["sherry"] # Pause immediately after Sherry's analysis
    )

app = create_mission_graph()