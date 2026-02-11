// frontend/src/App.tsx

import React, { useState } from 'react';
import MissionHub from './components/MissionHub';

const App: React.FC = () => {
  // 状態の初期化
  const [apiKey, setApiKey] = useState<string>(localStorage.getItem('GEMINI_API_KEY') || '');

  // キーの保存処理
  const handleSetApiKey = (val: string): void => {
    localStorage.setItem('GEMINI_API_KEY', val);
    setApiKey(val);
  };

  // キーの削除処理
  const handleDelete = (): void => {
    localStorage.removeItem('GEMINI_API_KEY');
    setApiKey('');
  };

  return (
    <div className="App">
      <MissionHub 
        apiKey={apiKey} 
        setApiKey={handleSetApiKey}
        onDelete={handleDelete}
      />
    </div>
  );
};

export default App;