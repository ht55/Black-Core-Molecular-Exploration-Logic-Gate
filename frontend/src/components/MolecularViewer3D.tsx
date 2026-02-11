// frontend/src/components/MolecularViewer3D.tsx

import React, { useEffect, useRef, useState } from 'react';
import * as $3Dmol from '3dmol';

interface MolecularViewer3DProps {
  smiles?: string;
  baseSmiles?: string;
  mode: '2d' | '3d'; 
}

declare global {
  interface Window {
    initRDKitModule?: () => Promise<any>;
    RDKit?: any;
  }
}

let rdkitInstance: any = null;

const MolecularViewer3D: React.FC<MolecularViewer3DProps> = ({ smiles, baseSmiles, mode }) => {
  const containerRef = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<any>(null);
  const [svgContent, setSvgContent] = useState<string>('');
  const [isRdkitReady, setIsRdkitReady] = useState(false);
  const [, forceUpdate] = useState({});
  const [transform, setTransform] = useState({ x: 0, y: 0, scale: 0.8 });
  const [isDragging, setIsDragging] = useState(false);
  const [lastMousePos, setLastMousePos] = useState({ x: 0, y: 0 });

  // 1. Initialize RDKit (既存のロジックを維持)
  useEffect(() => {
    if (rdkitInstance) {
      setIsRdkitReady(true);
      return;
    }
    const loadRDKit = async () => {
      if (window.initRDKitModule) {
        try {
          rdkitInstance = await window.initRDKitModule();
          setIsRdkitReady(true);
          forceUpdate({});
        } catch (e) {
          console.error("RDKit WASM init error:", e);
        }
      }
    };
    const timer = setInterval(() => {
      if (window.initRDKitModule) {
        clearInterval(timer);
        loadRDKit();
      }
    }, 200);
    return () => clearInterval(timer);
  }, []);

  // 2. 3D Visualization (既存のロジックを尊重しつつ、3Dモード時のみ実行)
  useEffect(() => {
    // 3Dモードでない場合や、依存値が足りない場合は実行しない
    if (mode !== '3d' || !containerRef.current || !rdkitInstance || !baseSmiles) return;

    const renderModel = () => {
      try {
        if (!viewerRef.current) {
          viewerRef.current = ($3Dmol as any).createViewer(containerRef.current, { 
            backgroundColor: 'black' 
          });
        }
        const v = viewerRef.current;
        v.clear();

        const bMol = rdkitInstance.get_mol(baseSmiles);

        if (bMol) {
          try {
            // エラーの出ないメソッド名（set_3d_coords等）を確認して実行
            if (typeof bMol.set_3d_coords === 'function') {
              bMol.set_3d_coords();
            } else if (typeof bMol.initialize_3d_coords === 'function') {
              bMol.initialize_3d_coords();
            }

            if (typeof bMol.mmff_has_mmff_params === 'function' && bMol.mmff_has_mmff_params()) {
              bMol.mmff_optimize();
            }

            // Prioritize V3000 & fallback
            const molData = bMol.get_v3000 ? bMol.get_v3000() : bMol.get_molblock();
            const format = bMol.get_v3000 ? "sdf" : "mol";

            v.addModel(molData, format);
            
            v.setStyle({}, { 
              stick: { color: 'spectrum', radius: 0.18 } 
            });

            v.zoomTo();
            v.render();

            setTimeout(() => {
              v.resize();
              v.zoomTo();
              v.render();
            }, 150);

          } catch (e) {
            console.error("3D_COORDS_OR_OPTIMIZATION_ERROR:", e);
            if (typeof bMol.set_2d_coords === 'function') bMol.set_2d_coords();
            v.addModel(bMol.get_molblock(), "mol");
            v.render();
          } finally {
            bMol.delete();
          }
          console.log("!!! 3D_V3000_MODEL_RENDERED !!!");
        }
      } catch (err) {
        console.error("3D_RENDER_ERROR:", err);
      }
    };

    requestAnimationFrame(renderModel);
  }, [baseSmiles, isRdkitReady, mode]);

  // 3. 2D Visualization 
  useEffect(() => {
    if (mode !== '2d' || !rdkitInstance || !smiles) return;
    try {
      const m = rdkitInstance.get_mol(smiles);
      if (m) {
        const svg = m.get_svg_with_highlights(JSON.stringify({
          width: 150, height: 150,
          bondLineWidth: 1.5,
          prepareMolsBeforeDrawing: true
        }));
        setSvgContent(svg);
        m.delete();
      }
    } catch (e) { console.error("2D_ERROR:", e); }
  }, [smiles, isRdkitReady, mode]);


  // --- Branching of the final rendering results ---
  if (mode === '2d') {
    return (
      <div 
        style={{ width: '100%', display: 'flex', justifyContent: 'center', alignItems: 'center' }}>
        {svgContent ? (
          <div 
            dangerouslySetInnerHTML={{ __html: svgContent }} 
            style={{ filter: 'invert(1) sepia(1) saturate(10) hue-rotate(90deg) brightness(1.2)' }} 
          />
        ) : (
          <div style={{ fontSize: '10px', color: '#005500', textAlign: 'center' }}>
            {!rdkitInstance ? "SYSTEM_OFFLINE\nINITIALIZING_WASM..." : "LOADING_GEOMETRY..."}
          </div>
        )}
      </div>
    );
  }

  // Return only the 3D-viz container
  return (
    <div ref={containerRef} style={{ width: '100%', height: '100%' }} />
  );
};

export default MolecularViewer3D;