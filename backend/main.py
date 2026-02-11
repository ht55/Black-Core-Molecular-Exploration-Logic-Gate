# backend/main.py

import os
import json
from fastapi import FastAPI, Query, Body, HTTPException
from fastapi.responses import StreamingResponse, FileResponse
from fastapi.middleware.cors import CORSMiddleware
from langchain_google_genai import ChatGoogleGenerativeAI
from graph_engine import app as chain
from engine import BlackPhysicsEngine
from dotenv import load_dotenv
from graph_engine import MASTER_ASSETS
from fastapi.staticfiles import StaticFiles

app = FastAPI()
load_dotenv()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

engine = BlackPhysicsEngine()

@app.post("/api/mission/stream/{mission_id}")
async def stream_mission(
    mission_id: str, 
    key: str = Query(None),
    payload: dict = Body(...)
):

    current_status = payload.get("status")
    lang = payload.get("language", "ja") 
    is_jp = (lang == "ja")

    if current_status == "FAILED":
        msg = "【API制限】APIリソースが尽きました。12:00 AM PT にリセット予定。明日またお試しください。" if is_jp else "[API LIMIT] API resources exhausted. Reset scheduled for 12:00 AM PT. Please try again tomorrow."
        raise HTTPException(status_code=429, detail=msg)

    api_key = key or os.getenv("GOOGLE_API_KEY")
    if not api_key or len(api_key) < 10:
        raise HTTPException(
            status_code=401, 
            detail="Unauthorized access. API Key is missing or invalid."
        )

    # LLM definition (Optimized for chemical structural modification)
    llm = ChatGoogleGenerativeAI(
        model="gemini-2.5-flash", 
        google_api_key=api_key,
        temperature=0.2, 
        max_output_tokens=5500,
        streaming=True
    )

    # payloadの中身を確定
    m_assets = MASTER_ASSETS.get(mission_id, {})
    scaffolds = m_assets.get("scaffolds", [])
    selected_id = payload.get("selected_scaffold_id", "")
    
    # --- base_scaffold_data の確定ロジック ---
    # 初回（iteration 0）またはデータが欠落している場合のみ、MASTER_ASSETSから引き当てる
    if payload.get("iteration", 0) == 0 or not payload.get("base_scaffold_data"):
        m_assets = MASTER_ASSETS.get(mission_id, {})
        scaffolds = m_assets.get("scaffolds", [])
        selected_id = payload.get("selected_scaffold_id", "")
        
        # IDが一致するものを引き当てる
        base_data = next((s for s in scaffolds if s["id"] == selected_id), {})
        payload["base_scaffold_data"] = base_data
    # iteration > 0 の場合は、フロントエンドから送られてきた base_scaffold_data をそのまま維持する

    # State Mapping for "Black-Hole Discovery" Architecture
    # フロントエンドからの新パラメータを優先的にマッピング
    state_input = {
        "mission_id": mission_id,
        "language": lang,
        "philosophy": payload.get("philosophy", "STANDARD").upper(),
        "scaffold_id": payload.get("selected_scaffold_id", ""),
        "base_scaffold_data": payload.get("base_scaffold_data", {}),
        "entropy_level": float(payload.get("entropy_level", 0.1)),
        
        "iteration": payload.get("iteration", 0),
        "status": "RUNNING", # ボタンが押されたらRUNNINGで再開
        "physics_fail_count": payload.get("physics_fail_count", 0),
        "current_smiles": payload.get("current_smiles", ""),
        "candidates": payload.get("candidates", []),
        "metrics_results": payload.get("metrics_results", []),
        "transformation_metrics": payload.get("transformation_metrics", {}),
        
        "scientific_report": payload.get("scientific_report", ""),
        "next_instruction": "",
        "dialogue": payload.get("dialogue", []),
        # "three_d_data": payload.get("three_d_data", {}),
        "distilled_lesson": payload.get("distilled_lesson", "")
    }

    # 内部ループの救済処置: physics_fail_count=3でNEXTが選択された場合、カウントを0に戻しphysics_nodeで即座にFAILEDになるのを防止。
    # if state_input["physics_fail_count"] >= 3:
        # state_input["physics_fail_count"] = 0

    async def event_generator():
        current_running_state = state_input.copy()
        config = {"configurable": {"llm": llm, "engine": engine, "thread_id": mission_id, "language": lang}}

        # デバッグ用
        if state_input["iteration"] > 0:
            print(f"DEBUG: [RE-STARTING FROM CHECKPOINT] Mission: {mission_id}, Iteration: {state_input['iteration']}")

        try:
            # iteration > 0 の場合の入力を整理
            initial_input = None if state_input["iteration"] > 0 else state_input

            # 🚀 
            yield f"data: {json.dumps({'type': 'node_start', 'node': 'vermouth'}, ensure_ascii=False)}\n\n"

            async for event in chain.astream(initial_input, config=config, stream_mode="updates"):
                print(f"DEBUG: Event received from chain: {event}") # for debagging
                for node_name, output in event.items():
                    if node_name.startswith("__"): continue

                    # 既存のStateを更新
                    current_running_state.update(output)

                    # ゾンビループ停止
                    # is_failed = current_running_state.get("status") == "FAILED"

                    # scientific_reportからシェリーのセリフに含める部分を抽出
                    if node_name == "sherry":
                        report = output.get("scientific_report", "")
                        metrics = output.get("metrics_results", [])
                        best_candidate = metrics[0] if metrics else {}
                        descriptors = best_candidate.get("descriptors", {})

                        # 2. Center Pane用の会話文を抽出
                        import re
                        display_text = ""
                        if report and report != "Analysis failed.":
                            # 最初の「\n\n」または「。 」で区切って、人間らしい結論だけを抽出
                            # data_separator_pattern = re.split(r'\n\n|(?<=。)\s', report)
                            paragraphs = report.split("\n\n")
                            display_text = paragraphs[0].strip()
                        
                        # 3. ステートの整合性と履歴の鎖を維持
                        if "dialogue" in output:
                            history = list(current_running_state.get("dialogue", []))
                            # new_responses = []
                            for msg in output["dialogue"]:
                                # シェリーの詳細見解がある場合、ここでのみ反映
                                if msg["agent"] == "Sherry":
                                    msg["text"] = display_text
                                
                                # 重複排除しつつ追加
                                if not any(h["agent"] == msg["agent"] and h["text"] == msg["text"] for h in history[-3:]):
                                    history.append(msg)
                            
                            # current_running_state と output の両方に最新の歴史を反映
                            current_running_state["dialogue"] = history
                            output["dialogue"] = history

                            # 万が一、dialogueにSherryがいなくてreportがある場合の救済
                            if not any(m["agent"] == "Sherry" for m in output["dialogue"]) and display_text:
                                output["dialogue"].insert(0, {"agent": "Sherry", "text": display_text})     

                        # 4. フロントエンドが直接参照できるキーを追加
                        # これにより、フロントエンド側で 'output.display_metrics' を見るだけで表示可能になる
                        output["display_metrics"] = descriptors
                        output["full_report"] = report
                        # output["current_sdf"] = best_candidate.get("sdf", "")

                        mut_val = best_candidate.get("mutation_score", 0) * 100
                        output["formatted_analysis"] = (
                            f"MUTATION_SCORE:{mut_val}% / MW:{descriptors.get('mw')} / PSA:{descriptors.get('psa')} / RINGS:{descriptors.get('rings')}"
                        )  

                    # デバッグ用にデータの整合性チェック
                    # print(f"DEBUG: [Sherry Node Trace] 3D Data Presence: {'three_d_data' in current_running_state}")
                    print(f"DEBUG: [Sherry Node Trace] Metrics Presence: {'metrics_results' in current_running_state}")
                    
                    response_payload = {
                        "type": "node_update",
                        "node": node_name,
                        "output": output,
                        "full_state": {**current_running_state}
                    }
                    yield f"data: {json.dumps(response_payload, ensure_ascii=False, default=str)}\n\n"
 
                    # 🚀
                    next_node = None
                    if current_running_state.get("status") == "FAILED":
                        print(f"DEBUG: Mission FAILED at {node_name}. Terminating stream immediately.")
                        break

                    if node_name == "vermouth":
                        if current_running_state.get("philosophy") == "SERENDIPITY":
                            next_node = "Mutate Engine" 
                        else:
                            next_node = "Physics Engine" 

                    elif node_name == "mutate":
                        next_node = "Physics Engine"

                    elif node_name == "physics":
                        next_node = "Sherry"
                        
                    if next_node:
                        yield f"data: {json.dumps({'type': 'node_start', 'node': next_node})}\n\n"

                # ブラウザのバッファをフラッシュさせるためのダミーを送信
                yield "data: {\"type\": \"end\"}\n\n"
                print("--- [STREAM FINISHED] Final state sent to frontend ---")

        except Exception as e:
            error_data = {"type": "error", "message": str(e)}
            yield f"data: {json.dumps(error_data)}\n\n"

    return StreamingResponse(event_generator(), media_type="text/event-stream")

# Mount the frontend build directory
build_path = os.path.join(os.getcwd(), "build")

if os.path.exists(build_path):
    # 1. Mount /static for JS, CSS, and hashed images (like blbg4.png)
    app.mount("/static", StaticFiles(directory=os.path.join(build_path, "static")), name="static")

    # 2. IMPORTANT: Mount the images directory explicitly if it exists in your build
    # This fixes the missing patterns in your SVG (imgNormal/imgHover)
    images_path = os.path.join(build_path, "images")
    if os.path.exists(images_path):
        app.mount("/images", StaticFiles(directory=images_path), name="images")

    @app.get("/{catchall:path}")
    async def read_index(catchall: str):
        # Serve any specific file from the build root (favicon, etc.)
        file_path = os.path.join(build_path, catchall)
        if os.path.isfile(file_path):
            return FileResponse(file_path)
        
        # Fallback to index.html for React SPA
        return FileResponse(os.path.join(build_path, "index.html"))