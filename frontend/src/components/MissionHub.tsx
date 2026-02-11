// frontend/src/components/MissionHub.tsx

import React, { useState } from 'react';
import './MissionHub.css';
import bgImage from './blbg4.png';
import MissionWindow from '../MissionWindow';

const container = document.getElementById('logo-container');

// --- Types & Interfaces ---
interface Mission {
  id: string;
  title: string;
  sub: string;
}

interface PathData {
  d: string;
  l1: string;
  l2: string;
}

// Appから渡されるデータの定義
interface MissionHubProps {
  apiKey: string;
  setApiKey: (val: string) => void;
  onDelete: () => void;
}

// model-viewerタグをTypeScriptに認識させるための宣言
declare global {
  namespace JSX {
    interface IntrinsicElements {
      'model-viewer': any;
    }
  }
}

const imgNormal = "/images/emeraldbg.png";
const imgHover = "/images/limebg.png";

const MissionHub: React.FC<MissionHubProps> = ({ apiKey, setApiKey, onDelete }) => {
  const [streamingText, setStreamingText] = useState("");
  const [isGenerating, setIsGenerating] = useState(false);

  const [isWindowOpen, setIsWindowOpen] = useState(false);
  const [selectedMission, setSelectedMission] = useState<Mission | null>(null);
  const [completedMissions, setCompletedMissions] = useState<string[]>([]);
  const [showKey, setShowKey] = useState<boolean>(false);

  const missions = [
    { id: 'm1', title: 'MISSION 001:', sub: 'STEALTH_BLACK' },
    { id: 'm2', title: 'MISSION 002:', sub: 'VOID_BLACK' },
    { id: 'm3', title: 'MISSION 003:', sub: 'THERMAL_BLACK' },
    { id: 'm4', title: 'MISSION 004:', sub: 'META_BLACK' }
  ];

  // パネルの形状と、2段テキスト用のパス（外側L1、内側L2）
  const paths = [
    {
      d: "M 310,25 A 275,275 0 0,1 575,290 L 435,290 A 135,135 0 0,0 310,165 Z",
      l1: "M 310,75 A 225,225 0 0,1 525,290",
      l2: "M 310,120 A 180,180 0 0,1 480,290"
    },
    { // MISSION 002 (右下) - 元を反転
      d: "M 575,310 A 275,275 0 0,1 310,575 L 310,435 A 135,135 0 0,0 435,310 Z",
      l1: "M 310,480 A 180,180 0 0,0 480,310",
      l2: "M 310,525 A 225,225 0 0,0 525,310"
    },
    { // MISSION 003 (左下) - 元を反転
      d: "M 290,575 A 275,275 0 0,1 25,310 L 165,310 A 135,135 0 0,0 290,435 Z",
      l1: "M 120,310 A 180,180 0 0,0 290,480",
      l2: "M 75,310 A 225,225 0 0,0 290,525"
    },
    {
      d: "M 25,290 A 275,275 0 0,1 290,25 L 290,165 A 135,135 0 0,0 165,290 Z",
      l1: "M 75,290 A 225,225 0 0,1 290,75",
      l2: "M 120,290 A 180,180 0 0,1 290,120"
    }
  ];

  // 入力時の処理  
  const handleSaveKey = () => {
    if (!apiKey) {
      alert("Please enter your Google API KEY.");
      return;
    }
    
    // Guidance alert in English
    alert("🔑 ACCESS KEY LOADED.\nProtocol ready. Key will be cleared when tab is closed.");
  };

  const handlePanelClick = (mission: Mission) => {
    if (!apiKey) {
      alert("Google API KEYを入力してください。解析プロトコルを開始できません。");
      return; // キーがない場合は開かないようにガード
    }
    setSelectedMission(mission);
    setIsWindowOpen(true);
    console.log(`[SYSTEM] Starting analysis for: ${mission.title}`);
  };

  return (
    <div
      className={`hub-root ${isWindowOpen ? 'is-dimmed' : ''}`}
      style={{
        backgroundImage: `url(${bgImage})`,
        backgroundSize: 'cover',
        backgroundPosition: '60% 30%',
        backgroundRepeat: 'no-repeat',
        position: 'relative' 
      }}
    >
      <div className="radial-layout-container">
        <svg 
          viewBox="0 0 620 620" 
          className="main-ui-svg" 
          preserveAspectRatio="xMidYMid meet"
          style={{ overflow: 'visible' }}
        >
          <defs>
            {/* 通常時の画像パターン */}
            <pattern id="pat-normal" patternUnits="userSpaceOnUse" width="300" height="300">
              <image href={imgNormal} x="0" y="0" width="300" height="300" preserveAspectRatio="none" />
            </pattern>

            {/* ホバー時の画像パターン */}
            <pattern id="pat-hover" patternUnits="userSpaceOnUse" width="300" height="300">
              <image href={imgHover} x="0" y="0" width="300" height="300" preserveAspectRatio="none" />
            </pattern>

            {/* 既存のフィルタ設定など */}
            {/* 浮かび上がらせるためのシャドウ/グロウ定義 */}
            <filter id="panel-glow" x="-20%" y="-20%" width="140%" height="140%">
              {/* 背景の暗いシャドウ */}
              <feDropShadow dx="0" dy="0" stdDeviation="5" floodColor="black" floodOpacity="0.8"/>
              {/* 外側に漏れる光（Glow） */}
              <feDropShadow dx="0" dy="0" stdDeviation="8" floodColor="emerald" floodOpacity="0.4"/>
            </filter>
          </defs>

          {missions.map((m, i) => (
            <g 
              key={m.id} 
              className="sector-group"
              // 解析を開始
              onClick={() => handlePanelClick(m)} 
              style={{ cursor: apiKey ? 'pointer' : 'not-allowed' }}
            >
              <path className="sector-arc-path" d={paths[i].d} />

              {/* 外側：MISSION: */}
              <path id={`path-t1-${i}`} d={paths[i].l1} fill="none" />
              <text className="m-text-main">
                <textPath href={`#path-t1-${i}`} startOffset="50%" textAnchor="middle">
                  {m.title}
                </textPath>
              </text>

              {/* 内側：MISSION詳細 */}
              <path id={`path-t2-${i}`} d={paths[i].l2} fill="none" />
              <text className="m-text-sub">
                <textPath href={`#path-t2-${i}`} startOffset="50%" textAnchor="middle">
                  {m.sub}
                </textPath>
              </text>
            </g>
          ))}
        </svg>

        <div className="center-command-unit">
          <div className="unit-display" style={{ textAlign: 'center' }}>
            <div className="unit-items-wrapper">
              {/* 3Dロゴ：model-viewerを使用したメインビジュアル */}
              <div className="logo-3d-container">
                <model-viewer
                  src="/images/titlelogo.gltf"
                  alt="Title Logo"
                  auto-rotate
                  camera-controls
                  disable-zoom
                  style={{ width: '100%', height: '100%', '--poster-color': 'transparent' }}
                  shadow-intensity="1"
                  environment-image="neutral"
                  exposure="1"
                />
              </div>

              {/* 一言メモ：Black Molecular Research */}
              <h4 className="unit-memo orbitron">
                  Molecular Exploration<br />
                  Logic Gate
              </h4>

              {/* 認証エリア：ACCESS KEY入力欄 */}
              <div className="unit-auth-area">
                <input 
                  type={showKey ? "text" : "password"}
                  value={apiKey}
                  onChange={(e) => {
                    setApiKey(e.target.value);
                  }}
                  placeholder="🔒 Google API KEY"
                />
                <div className="unit-controls">
                  <span onClick={() => setShowKey(!showKey)}>{showKey ? 'HIDE' : 'SHOW'}</span>
                  {/* ERASE時もメモリから消すだけにする */}
                  <span onClick={() => { setApiKey(''); }}>ERASE</span>
            
                  <a href="https://aistudio.google.com/app/apikey" target="_blank" rel="noreferrer" className="unit-auth-btn">GET</a>
                </div>
              </div>

              {/* 実行ボタン：状態によって色が変わる */}
              <button 
                className={`unit-execute-btn ${apiKey ? 'active' : ''}`}
                onClick={handleSaveKey}
              >
                {apiKey ? 'LOAD KEY' : 'LOCKED'}
              </button>
            </div>
          </div>
        </div>
      </div>

      {/* ウィンドウコンポーネント: selectedMissionのnullチェックを追加し、ロジックはそのまま維持 */}
      {isWindowOpen && selectedMission && (
        <MissionWindow
          key={selectedMission.id} 
          isOpen={isWindowOpen}
          missionData={selectedMission}
          apiKey={apiKey} 
          onClose={(wasSuccess?: any) => {
                      if (wasSuccess === true) {
                        setCompletedMissions(prev => {
                          // includesがエラーになる場合は (prev.indexOf(selectedMission.id) !== -1) にしてください
                          if (prev.includes(selectedMission.id)) return prev;
                          return [...prev, selectedMission.id];
                        });
                      }
                      setIsWindowOpen(false);
                    }}
                  />
                )}
              </div>
            );
          };

export default MissionHub;