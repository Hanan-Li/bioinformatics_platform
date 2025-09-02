'use client';

import { useState, useCallback, useRef } from 'react';
import Navigation from '../components/Navigation';

type NodeType = {
  id: string;
  type: string;
  position: { x: number; y: number };
  data: { label: string; description: string };
};

type EdgeType = {
  id: string;
  source: string;
  target: string;
};

const operationTypes = [
  { type: 'input', label: 'Data Input', description: 'Input data files', color: 'bg-blue-500' },
  { type: 'quality_control', label: 'Quality Control', description: 'QC analysis', color: 'bg-green-500' },
  { type: 'alignment', label: 'Alignment', description: 'Sequence alignment', color: 'bg-yellow-500' },
  { type: 'variant_calling', label: 'Variant Calling', description: 'Call genetic variants', color: 'bg-purple-500' },
  { type: 'annotation', label: 'Annotation', description: 'Annotate variants', color: 'bg-red-500' },
  { type: 'output', label: 'Output', description: 'Export results', color: 'bg-gray-500' },
];

export default function CreateWorkflow() {
  const [nodes, setNodes] = useState<NodeType[]>([]);
  const [edges, setEdges] = useState<EdgeType[]>([]);
  const [selectedNode, setSelectedNode] = useState<string | null>(null);
  const [draggedNode, setDraggedNode] = useState<string | null>(null);
  const canvasRef = useRef<HTMLDivElement>(null);
  const [isConnecting, setIsConnecting] = useState<string | null>(null);

  const onDragStart = (event: React.DragEvent, nodeType: string) => {
    event.dataTransfer.setData('application/reactflow', nodeType);
    event.dataTransfer.effectAllowed = 'move';
  };

  const onDrop = useCallback(
    (event: React.DragEvent) => {
      event.preventDefault();

      if (!canvasRef.current) return;

      const reactFlowBounds = canvasRef.current.getBoundingClientRect();
      const type = event.dataTransfer.getData('application/reactflow');

      if (typeof type === 'undefined' || !type) {
        return;
      }

      const position = {
        x: event.clientX - reactFlowBounds.left,
        y: event.clientY - reactFlowBounds.top,
      };

      const operationType = operationTypes.find(op => op.type === type);
      if (!operationType) return;

      const newNode: NodeType = {
        id: `${type}_${Date.now()}`,
        type,
        position,
        data: {
          label: operationType.label,
          description: operationType.description,
        },
      };

      setNodes((nds) => nds.concat(newNode));
    },
    []
  );

  const onDragOver = (event: React.DragEvent) => {
    event.preventDefault();
    event.dataTransfer.dropEffect = 'move';
  };

  const onNodeClick = (nodeId: string) => {
    if (isConnecting) {
      if (isConnecting !== nodeId) {
        const newEdge: EdgeType = {
          id: `e${isConnecting}-${nodeId}`,
          source: isConnecting,
          target: nodeId,
        };
        setEdges((eds) => eds.concat(newEdge));
      }
      setIsConnecting(null);
    } else {
      setSelectedNode(nodeId === selectedNode ? null : nodeId);
    }
  };

  const startConnection = (nodeId: string) => {
    setIsConnecting(nodeId);
    setSelectedNode(null);
  };

  const deleteNode = (nodeId: string) => {
    setNodes((nds) => nds.filter((node) => node.id !== nodeId));
    setEdges((eds) => eds.filter((edge) => edge.source !== nodeId && edge.target !== nodeId));
    setSelectedNode(null);
  };

  const deleteEdge = (edgeId: string) => {
    setEdges((eds) => eds.filter((edge) => edge.id !== edgeId));
  };

  const saveWorkflow = () => {
    const workflow = {
      nodes,
      edges,
      created_at: new Date().toISOString(),
    };
    console.log('Saving workflow:', workflow);
    alert('Workflow saved! (Check console for details)');
  };

  return (
    <div className="min-h-screen bg-gray-50">
      <Navigation />
      
      <div className="flex h-screen pt-16">
        <div className="w-64 bg-white border-r border-gray-200 p-4">
          <h2 className="text-lg font-semibold mb-4">Operations</h2>
          <div className="space-y-2">
            {operationTypes.map((op) => (
              <div
                key={op.type}
                className={`${op.color} text-white p-3 rounded-lg cursor-grab hover:opacity-80 transition-opacity`}
                draggable
                onDragStart={(event) => onDragStart(event, op.type)}
              >
                <div className="font-medium text-sm">{op.label}</div>
                <div className="text-xs opacity-90">{op.description}</div>
              </div>
            ))}
          </div>

          <div className="mt-8">
            <h3 className="text-md font-semibold mb-2">Actions</h3>
            <button
              onClick={saveWorkflow}
              className="w-full bg-blue-600 text-white p-2 rounded-lg hover:bg-blue-700 transition-colors mb-2"
            >
              Save Workflow
            </button>
            <button
              onClick={() => {
                setNodes([]);
                setEdges([]);
                setSelectedNode(null);
              }}
              className="w-full bg-red-600 text-white p-2 rounded-lg hover:bg-red-700 transition-colors"
            >
              Clear Canvas
            </button>
          </div>

          {selectedNode && (
            <div className="mt-6 p-3 bg-gray-100 rounded-lg">
              <h3 className="font-medium mb-2">Selected Node</h3>
              <p className="text-sm text-gray-600 mb-2">ID: {selectedNode}</p>
              <button
                onClick={() => startConnection(selectedNode)}
                className="w-full bg-green-600 text-white p-1 rounded text-sm hover:bg-green-700 mb-2"
              >
                Connect from here
              </button>
              <button
                onClick={() => deleteNode(selectedNode)}
                className="w-full bg-red-600 text-white p-1 rounded text-sm hover:bg-red-700"
              >
                Delete Node
              </button>
            </div>
          )}

          {isConnecting && (
            <div className="mt-4 p-3 bg-yellow-100 rounded-lg">
              <p className="text-sm">Click another node to connect</p>
              <button
                onClick={() => setIsConnecting(null)}
                className="mt-2 text-sm text-red-600 hover:text-red-800"
              >
                Cancel
              </button>
            </div>
          )}
        </div>

        <div className="flex-1 relative">
          <div
            ref={canvasRef}
            className="w-full h-full bg-gray-100"
            onDrop={onDrop}
            onDragOver={onDragOver}
            style={{
              backgroundImage: 'radial-gradient(circle, #d1d5db 1px, transparent 1px)',
              backgroundSize: '20px 20px'
            }}
          >
            <svg className="absolute inset-0 w-full h-full pointer-events-none">
              {edges.map((edge) => {
                const sourceNode = nodes.find(n => n.id === edge.source);
                const targetNode = nodes.find(n => n.id === edge.target);
                
                if (!sourceNode || !targetNode) return null;
                
                return (
                  <line
                    key={edge.id}
                    x1={sourceNode.position.x + 50}
                    y1={sourceNode.position.y + 25}
                    x2={targetNode.position.x + 50}
                    y2={targetNode.position.y + 25}
                    stroke="#374151"
                    strokeWidth="2"
                    markerEnd="url(#arrowhead)"
                  />
                );
              })}
              <defs>
                <marker
                  id="arrowhead"
                  markerWidth="10"
                  markerHeight="7"
                  refX="9"
                  refY="3.5"
                  orient="auto"
                >
                  <polygon
                    points="0 0, 10 3.5, 0 7"
                    fill="#374151"
                  />
                </marker>
              </defs>
            </svg>

            {nodes.map((node) => {
              const operationType = operationTypes.find(op => op.type === node.type);
              return (
                <div
                  key={node.id}
                  className={`absolute w-24 h-16 ${operationType?.color || 'bg-gray-500'} text-white rounded-lg cursor-pointer shadow-md hover:shadow-lg transition-shadow ${
                    selectedNode === node.id ? 'ring-2 ring-blue-400' : ''
                  } ${isConnecting === node.id ? 'ring-2 ring-yellow-400' : ''}`}
                  style={{
                    left: node.position.x,
                    top: node.position.y,
                  }}
                  onClick={() => onNodeClick(node.id)}
                >
                  <div className="p-2 text-xs text-center">
                    <div className="font-medium truncate">{node.data.label}</div>
                  </div>
                  {edges.some(edge => edge.source === node.id || edge.target === node.id) && (
                    <div className="absolute -top-1 -right-1 w-3 h-3 bg-green-400 rounded-full"></div>
                  )}
                </div>
              );
            })}

            {nodes.length === 0 && (
              <div className="absolute inset-0 flex items-center justify-center">
                <div className="text-center text-gray-500">
                  <svg className="mx-auto h-12 w-12 text-gray-400" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M7 16a4 4 0 01-.88-7.903A5 5 0 1115.9 6L16 6a5 5 0 011 9.9M15 13l-3-3m0 0l-3 3m3-3v12" />
                  </svg>
                  <h3 className="mt-2 text-sm font-medium">No workflow nodes</h3>
                  <p className="mt-1 text-sm">Drag operations from the sidebar to get started.</p>
                </div>
              </div>
            )}
          </div>
        </div>
      </div>
    </div>
  );
}