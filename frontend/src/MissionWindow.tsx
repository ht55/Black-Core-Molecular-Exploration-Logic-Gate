// frontend/src/MissionWindow.tsx

import React, { useState, useEffect, useRef } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import MolecularViewer3D from './components/MolecularViewer3D';

// 追加する型定義
interface Scaffold {
  id: string;
  name: string;
  desc: string;
  smiles: string;
}

interface MissionAsset {
  label: string;
  scaffolds: Scaffold[];
}

// ミッションごとのScaffold定義（バックエンドのSCAFFOLD_LIBRARYと同期）
const MISSION_ASSETS: Record<string, MissionAsset> = {
  "m1": { 
    label: "Stealth Black", 
    scaffolds: [
      { id: "s1_1", name: "Azo-Linker", desc: "NIR Absorption Core", smiles: "C1=CC=C(C=C1)N=NC2=CC=CC=C2" },
      { id: "s1_2", name: "Phthalocyanine", desc: "Robust Photon Trap", smiles: "c1ccc2c(c1)c3nc4nc(nc5nc(nc6nc(n3)c2c6)c7ccccc75)c8ccccc84" }
    ]
  },
  "m2": { 
    label: "Void Black", 
    scaffolds: [
      { id: "s2_1", name: "Perylene-Core", desc: "Aromatic Labyrinth", smiles: "c1cc2cccc3c2c(c1)ccc4cccc3c4" },
      { id: "s2_2", name: "Coronene", desc: "Zero-Reflectance Target", smiles: "c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67" }
    ]
  },
  "m3": { 
    label: "Thermal Black", 
    scaffolds: [
      { id: "s3_1", name: "Triphenylene", desc: "Discotic Liquid Crystal", smiles: "c1ccc2c(c1)c3ccc4ccccc4c3c5ccccc25" },
      { id: "s3_2", name: "Anthracene-Dimer", desc: "PCM Base Structure", smiles: "c1ccc2c(c1)cc3ccc4ccccc4c3c2" }
    ]
  },
  "m4": { 
    label: "Meta Black", 
    scaffolds: [
      { id: "s4_1", name: "Helicene", desc: "Wavefront Manipulation", smiles: "c1ccc2c(c1)ccc3c2ccc4c3ccc5c4ccc6ccccc65" },
      { id: "s4_2", name: "D-A-Quaterthiophene", desc: "Electronic Camouflage", smiles: "s1cccc1c2sccc2c3sccc3c4sccc4" }
    ]
  },
};

interface LogItem {
  id: string;
  agent: string;
  text: string;
  critical: boolean;
}

interface AnalysisState {
  val: string | number;
  unit: string;
}

interface Props {
  isOpen: boolean;
  missionData: { id: string; sub: string };
  onClose: () => void;
  apiKey: string;
}

const MissionWindow: React.FC<Props> = ({ isOpen, missionData, onClose, apiKey }) => {
  const [logs, setLogs] = useState<LogItem[]>([]);
  const [status, setStatus] = useState<string>('IDLE');
  const [currentSmiles, setCurrentSmiles] = useState<string>("");
  const [lang, setLang] = useState<'ja' | 'en'>('en');
  const [analysis, setAnalysis] = useState<AnalysisState>({ val: "--", unit: "%" });
  // const [sdfData, setSdfData] = useState<string | null>(null);
  
  const [philosophy, setPhilosophy] = useState<string>('STANDARD');
  const [selectedScaffold, setSelectedScaffold] = useState<string>("");
  const [entropy, setEntropy] = useState<number>(0.1); 

  const scrollRef = useRef<HTMLDivElement>(null);
  const isRunningRef = useRef<boolean>(false);
  const abortControllerRef = useRef<AbortController | null>(null);
  const [fullState, setFullState] = useState<any>(null);
  const shownMessageIdsRef = useRef<Set<string>>(new Set());
  const [activeNode, setActiveNode] = useState<string | null>(null);

  useEffect(() => {
    if (missionData && MISSION_ASSETS[missionData.id]) {
      setSelectedScaffold(MISSION_ASSETS[missionData.id].scaffolds[0].id);
    }
  }, [missionData]);

  useEffect(() => {
    if (scrollRef.current) scrollRef.current.scrollTop = scrollRef.current.scrollHeight;
  }, [logs]);

  const typeOutLog = async (agent: string, text: string, isUpdate: boolean) => {
    const newLog = { 
      id: `${agent}-${Date.now()}-${Math.random()}`, 
      agent, 
      text,
      critical: isUpdate 
    };
    setLogs(prev => [...prev, newLog]);

    if (scrollRef.current) {
      window.requestAnimationFrame(() => {
        if (scrollRef.current) {
          scrollRef.current.scrollTop = scrollRef.current.scrollHeight;
        }
      });
    }
  };

  const runProtocol = async () => {
    console.log("LOG_4: runProtocol entered. isRunningRef.current is:", isRunningRef.current); //デバッグ用
    if (isRunningRef.current) {
      console.log("LOG_5: runProtocol ABORTED by isRunningRef guard");
      return;
    }

    abortControllerRef.current = new AbortController();
    isRunningRef.current = true;
    setStatus('RUNNING');
    
    if (!fullState) setLogs([]); 

    const payload = {
      ...(fullState || { iteration: 0, dialogue: [] }),
      mission_id: missionData.id,
      language: lang,
      // 最初だけ外側のphilosophy/entropyを使い、2回目以降（iteration > 0）は、引き継いだfullStateの値を尊重する
      philosophy: fullState?.iteration > 0 ? fullState.philosophy : philosophy,
      entropy_level: fullState?.iteration > 0 ? fullState.entropy_level : entropy,
      selected_scaffold_id: selectedScaffold,
    };

    try {
      console.log("SENDING_PAYLOAD:", payload); // for debagging
      const response = await fetch(`/api/mission/stream/${missionData.id}?key=${apiKey}`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(payload),
        signal: abortControllerRef.current.signal
      });

      // 429 Error response handling (API limit reached)
      if (!response.ok) {
        const errorData = await response.json().catch(() => ({}));
        alert(`MISSION FAILED\n\n${errorData.detail}`);
        
        isRunningRef.current = false;
        setStatus('IDLE');
        return; 
      }
      
      if (!response.body) return;
      const reader = response.body.getReader();
      const decoder = new TextDecoder();
      let buffer = "";

      while (true) {
        const { value, done } = await reader.read();
        if (done) break;
        buffer += decoder.decode(value, { stream: true });
        const lines = buffer.split("\n\n"); 
        buffer = lines.pop() || "";

        for (const line of lines) {
          const trimmedLine = line.trim();
          if (!trimmedLine.startsWith("data: ")) continue;

          try {
            const data = JSON.parse(trimmedLine.replace("data: ", ""));
            console.log("RECEIVED_DATA:", data); // for debagging

            // 🚀 Status logic
            if (data.type === "node_start") {
                setActiveNode(data.node);

                const nodeLabel = data.node.toUpperCase();
                typeOutLog("SYSTEM", `>> AWAITING [${nodeLabel}]...`, false);
                continue;
            }
            if (data.node) {
                setActiveNode(data.node);
            }

            const node = data.node;
            const output = data.output || {};
            console.log("RECEIVED_DATA:", data); // for debagging
            
            const currentFullState = data.full_state;

            if (currentFullState) {
              const updatedFullState = node === "sherry" 
                ? { ...currentFullState, show_user_choice: true } 
                : currentFullState;
              setFullState(updatedFullState);

              if (updatedFullState.current_smiles) {
                setCurrentSmiles(updatedFullState.current_smiles);
              }
              /*　if (updatedFullState.three_d_data) {
                const s_val = Object.values(updatedFullState.three_d_data)[0] as string;
                if (s_val) setSdfData(s_val);
              } */

              if (updatedFullState.status === 'FAILED') {
                setStatus('IDLE');
                break; 
              }
            }

            // From all nodes
            // 1. Design Node (Vermouth)
            if (data.node === "vermouth") {
              for (const msg of (output.dialogue || [])) {
                const msgKey = msg.id || `${msg.agent}-${msg.text}`;
                if (!shownMessageIdsRef.current.has(msgKey)) {
                  shownMessageIdsRef.current.add(msgKey);
                  await typeOutLog(msg.agent, msg.text, false);
                }
              }
            }
            
            // 2. Physics Node (物理エンジン)
            if (data.node === "physics engine") {
              // const sdfs = output.three_d_data || {};
              // const targetSdf = Object.values(sdfs)[0] as string;
              //if (targetSdf) {
                // setSdfData(targetSdf);
                //console.log("3D Data (SDF) Loaded");
              //}

              const results = output.metrics_results || [];
              if (results[0]) {
                setCurrentSmiles(results[0].smiles);

                const score = results[0].mutation_score;
                if (typeof score === 'number' && !isNaN(score)) {
                  setAnalysis({ val: (score * 100).toFixed(1), unit: "%" });
                }
              }

              for (const m of (output.dialogue || [])) {
                const key = m.id || `${m.agent}-${m.text}`;
                if (!shownMessageIdsRef.current.has(key)) {
                  shownMessageIdsRef.current.add(key);
                  await typeOutLog(m.agent, m.text, true);
                }
              }
            }
            
            // 3. Analysis Node (Sherry)
            if (data.node === "sherry") {
              console.log("Sherry Node Data:", data);

              // Get SMILES from output or full_state
              const s_smiles = output.current_smiles || currentFullState?.current_smiles;
              if (s_smiles) {
                setCurrentSmiles(s_smiles);

                // if (output.scientific_report) {
                // await typeOutLog("Sherry", `【SCIENTIFIC REPORT】\n${output.scientific_report}`, true);}
              }

              // sherryノード通過時に変異率表示が消えないよう再セット
              const s_results = output.metrics_results || currentFullState?.metrics_results || [];

              // ログに記録されている metrics_results[0].sdf を優先
              // const targetSdf = s_results[0]?.sdf || (output.three_d_data ? Object.values(output.three_d_data)[0] : null);
              // if (targetSdf && typeof targetSdf === 'string') {
                // setSdfData(targetSdf.trim()); // 前後の改行をトリムして確実にパースさせる
                // console.log("3D Data (SDF) Successfully set in Sherry Node");
              // }

              if (s_results[0] && s_results[0].mutation_score !== undefined) {
                setAnalysis({ val: (s_results[0].mutation_score * 100).toFixed(1), unit: "%" });
              }

              if (output.formatted_analysis) {
                setAnalysis({ 
                  val: output.formatted_analysis, 
                  unit: "" 
                });
              }  

              // ログ出力（hallucination blocker）
              for (const msg of output.dialogue || []) {
                const msgKey = `${msg.agent}-${msg.text.trim()}`;
                
                if (msg.text.trim() !== "" && !shownMessageIdsRef.current.has(msgKey)) {
                    shownMessageIdsRef.current.add(msgKey);
                    await typeOutLog(msg.agent, msg.text, false);
                    
                }
              }
            }
          } catch (e) { console.error("JSON_PARSE_ERROR", e); }
        }
      }
    } catch (err: any) {
      if (err.name !== 'AbortError') {
        console.error("STREAM_ERROR:", err);
        const msg = err.response?.data?.detail || err.message;
        alert(`MISSION FAILED\n\n${msg}`);

        setStatus('IDLE');
      }
    } finally {
      isRunningRef.current = false;
      setStatus('IDLE');
    }
  };

  const currentAsset = MISSION_ASSETS[missionData.id];

  // Next Iteration
  const handleNextIteration = async () => {
    console.log("LOG_1: handleNextIteration CLICKED"); //デバッグ用

    isRunningRef.current = false;
    setStatus('IDLE');

    console.log("LOG_2: isRunningRef set to false, status set to IDLE"); //デバッグ用

    setTimeout(() => {
      console.log("LOG_3: setTimeout callback executing. calling runProtocol...");
      runProtocol();
    }, 50);
  };

  // handleArchive (SAVE RESULT)
  const handleArchive = () => {
    const content = `[BLACK_ORGANIZATION_MISSION_REPORT]\n` +
                    `DATE: ${new Date().toLocaleString()}\n` +
                    `MISSION_ID: ${missionData.id}\n` +
                    `SCAFFOLD: ${selectedScaffold}\n` +
                    `PHILOSOPHY: ${philosophy}\n` +
                    `ENTROPY: ${(entropy * 100).toFixed(0)}%\n` +
                    `ITERATION: ${fullState?.iteration || 0}\n` +
                    `------------------------------------------\n` +
                    `RESULT_SMILES: ${currentSmiles || "NONE (Physics Violation)"}\n` +
                    `ANALYSIS_VALUE: ${analysis.val} ${analysis.unit}\n` +
                    `------------------------------------------\n` +
                    `SCIENTIFIC_REPORT:\n${fullState?.scientific_report || "NO_REPORT_AVAILABLE"}\n`;

    const blob = new Blob([content], { type: 'text/plain' });
    const url = URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href = url;
    link.download = `result_${missionData.id}_it${fullState?.iteration || 0}.txt`;
    link.click();
    URL.revokeObjectURL(url);

    setStatus('SUCCESS');
    onClose();
  };

return (
  <AnimatePresence>
    {isOpen && (
      <div className="fixed inset-0 flex items-center justify-center z-[9999] p-4 font-mono">
        <motion.div initial={{ opacity: 0 }} animate={{ opacity: 1 }} className="absolute inset-0 bg-black/90" onClick={onClose} />
        <motion.div className="glass-modal relative w-full h-full max-w-7xl border-2 border-green-500 flex flex-col overflow-hidden bg-black text-green-500">
          
          {/* mission-layout-rootと3つのpaneクラス */}
          <div className="mission-layout-root">
            
            {/* LEFT PANE: Settings & Actions */}
            <div className="pane-left">
              <div className="control-group">
                <div className="section-label">TARGET_MISSION:</div>
                <div className="mission-id">{missionData.id.replace('m', '').padStart(3, '0')}: {missionData.sub}</div>
              </div>

              {/* Base Scaffold */}
              <div className="control-group">
                <div className="section-label">BASE_SCAFFOLD:</div>
                <div className="parameter-note">Determines the starting core structure and mutation sites.</div>
                <select 
                  value={selectedScaffold} 
                  onChange={(e) => setSelectedScaffold(e.target.value)}
                  className="scafold-selection"
                >
                  {currentAsset?.scaffolds.map((s: any) => (
                    <option key={s.id} value={s.id}>{s.name}</option>
                  ))}
                </select>
              </div>

              {/* Philosophy */}
              <div className="control-group">
                <div className="section-label">EXPLORATION_PHILOSOPHY:</div>
                <div className="parameter-note">Defines the mutation logic: Conservative vs. Radical skeletal change.</div>
                <div className="flex flex-col gap-2"> 
                  {['Standard', 'Taboo', 'Serendipity'].map(p => (
                    <button 
                      key={p} 
                      onClick={() => setPhilosophy(p.toUpperCase())}
                      className={`btn-pill-list ${philosophy === p.toUpperCase() ? 'active' : ''}`}
                    >
                      <span className="arrow-icon">➤</span> 
                      {p}
                    </button>
                  ))}
                </div>
              </div>

              {/* Entropy Slider-bar */}
              <div className="control-group">
                <div className="section-label">ENTROPY_LEVEL: {(entropy * 100).toFixed(0)}%</div>
                <div className="parameter-note">Controls the intensity and quantity of structural modifications.</div>
                <input 
                  type="range" min="0" max="1" step="0.01"
                  value={entropy} 
                  onChange={(e) => setEntropy(parseFloat(e.target.value))} 
                  className="custom-slider" />
              </div>

              {/* Invoke Protocol button */}
              <button onClick={runProtocol} className="btn-execute active px-10 w-fit">INVOKE_PROTOCOL</button>

              {/* Choice after Sherry's analysis */}
              {fullState?.show_user_choice && (
                <div className="flex gap-2 animate-fade-in mt-4">
                  <button 
                    onClick={handleNextIteration} 
                    className="btn-next"
                    disabled={status === 'RUNNING'} 
                  >
                    NEXT<br />ITERATION
                  </button>
                  <button 
                    onClick={handleArchive} 
                    className="btn-save"
                    disabled={status === 'RUNNING'} 
                  >
                    SAVE<br />RESULT
                  </button>
                </div>
              )}  

              {/* Abort Mission button */}      
              <div className="left-bottom-actions">
                <button onClick={onClose} className="btn-secondary">ABORT_MISSION</button>
              </div>
            </div>

            {/* CENTER PANE */}
            <div ref={scrollRef} className="pane-center">
              {/* Language button: Upper left corner */}
              <div className="lang-selector-container">
                <div className="lang-button-group">
                  <button 
                    className={`lang-btn ${lang === 'en' ? 'active' : ''}`}
                    onClick={() => setLang('en')}
                  >
                    EN
                  </button>
                  <div className="lang-btn-divider" />
                  <button 
                    className={`lang-btn ${lang === 'ja' ? 'active' : ''}`}
                    onClick={() => setLang('ja')}
                  >
                    JA
                  </button>
                </div>
              </div>

              {/* Agents dialogue: ログ表示ロジック */}
              {(!fullState || !fullState.dialogue || fullState.dialogue.length === 0) && (
                <div className="section-label animate-pulse">AWAITING_INPUT...</div>
              )}
              
              {logs.map((log, index) => (
                <div key={log.id || index} className="log-entry">
                  <span className="agent-name" data-agent={log.agent}>[{log.agent}]</span>
                  <span className="message-text">{log.text}</span>
                </div>
              ))}
            </div>

            {/* RIGHT PANE: Realtime analysis & visualizations */}
            <div className="pane-right">
              <div className="section-label">REALTIME_ANALYSIS:</div>
              <div className="metric-num">
                {analysis.val} 
                <small className="text-sm ml-1">{analysis.unit}</small>
              </div>

              {/* NEW_MOLECULE (2D) section */}
              <div className="section-label">NEW_MOLECULE:</div>
              <div className="smiles-structure">
                {currentSmiles || "AWAITING_DATA..."}
              </div>
              {/* 2D-viz container */}
              <div className="viz-container is-2d">
                <MolecularViewer3D mode="2d" smiles={currentSmiles} />
              </div>

              {/* BASE_MOLECULE (3D) section */}
              <div className="section-label">BASE_MOLECULE:</div>
              <div className="basemol-note">Interactive 3D view.</div>
              <div className="viz-container is-3d">
                {(() => {
                  const targetScaffold = currentAsset?.scaffolds?.find((s: any) => s.id === selectedScaffold);
                  const activeBaseSmiles = targetScaffold?.smiles || "";
                  return <MolecularViewer3D mode="3d" baseSmiles={activeBaseSmiles} />;
                })()}
              </div>
            </div>
          </div>
        </motion.div>
      </div>
    )}
  </AnimatePresence>
  );
};

export default MissionWindow;