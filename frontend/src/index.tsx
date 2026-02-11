// frontend/src/index.tsx

import React from 'react';
import ReactDOM from 'react-dom/client';
import './components/MissionHub.css';
import App from './App';

// getElementByIdがnullを返す可能性を考慮し、!を付ける
const rootElement = document.getElementById('root');

if (!rootElement) throw new Error('Failed to find the root element');

const root = ReactDOM.createRoot(rootElement);
root.render(
  <React.StrictMode>
    <App />
  </React.StrictMode>
);