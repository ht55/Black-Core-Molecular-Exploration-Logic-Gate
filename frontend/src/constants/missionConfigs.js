// frontend/src/constants/missionConfigs.js

export const MISSION_CONFIGS = {
  m1: {
    id: "m1",
    label: "NIR_Reflectance",
    displayName: "近赤外線反射率",
    unit: "%",
    min: 0.01,   // 1% (究極のステルス)
    max: 0.10,   // 10% (ベースライン)
    step: 0.01,
    default: 0.05,
    autonomous_base: 0.05, // アーカイブの到達点
    description: "800nm以上の赤外線透過を維持しつつ、反射を極限まで抑えろ。"
  },
  m2: {
    id: "m2",
    label: "Visible_Reflectance",
    displayName: "可視光反射率",
    unit: "%",
    min: 0.001,  // 0.001% (狂気の暗黒)
    max: 0.05,   // 5% (一般的な黒)
    step: 0.001,
    default: 0.03,
    autonomous_base: 0.03, // Geminiでの到達点
    description: "光を閉じ込める迷宮を構築し、肉眼による認識を無効化せよ。"
  },
  m3: {
    id: "m3",
    label: "Emissivity",
    displayName: "放射率",
    unit: "ε",
    min: 0.01,   // 低放射 (Low-E)
    max: 1.00,   // 高放射 (背景同調用)
    step: 0.01,
    default: 0.95,
    autonomous_base: 0.95, // ログの目標値
    description: "周囲の背景放射率に同調し、サーマルカメラの目を欺け。"
  },
  m4: {
    id: "m4",
    label: "Refractive_Index",
    displayName: "屈折率",
    unit: "n",
    min: -2.0,   // 強力な負の屈折
    max: 1.0,    // 通常の屈折率(空気=1)
    step: 0.1,
    default: -1.0,
    autonomous_base: -1.0, // 阿笠博士の到達点
    description: "物理法則を歪め、光を迂回させるメタ構造を設計せよ。"
  }
};