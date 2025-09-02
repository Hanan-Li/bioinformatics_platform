'use client';

import { useState } from 'react';
import Navigation from '../components/Navigation';

type WorkflowItem = {
  id: string;
  name: string;
  status: 'completed' | 'running' | 'failed' | 'pending';
  created_at: string;
  nodes_count: number;
  description: string;
};

const mockWorkflows: WorkflowItem[] = [
  {
    id: '1',
    name: 'RNA-seq Analysis Pipeline',
    status: 'completed',
    created_at: '2024-01-15T10:30:00Z',
    nodes_count: 6,
    description: 'Complete RNA sequencing analysis from raw reads to differential expression'
  },
  {
    id: '2',
    name: 'Variant Calling Workflow',
    status: 'running',
    created_at: '2024-01-14T14:22:00Z',
    nodes_count: 4,
    description: 'Identify genetic variants from whole genome sequencing data'
  },
  {
    id: '3',
    name: 'ChIP-seq Analysis',
    status: 'failed',
    created_at: '2024-01-13T09:15:00Z',
    nodes_count: 8,
    description: 'Chromatin immunoprecipitation sequencing analysis pipeline'
  },
  {
    id: '4',
    name: 'Metagenomics Classification',
    status: 'pending',
    created_at: '2024-01-12T16:45:00Z',
    nodes_count: 5,
    description: 'Taxonomic classification of metagenomic samples'
  }
];

export default function ViewWorkflow() {
  const [workflows] = useState<WorkflowItem[]>(mockWorkflows);
  const [filter, setFilter] = useState<'all' | 'completed' | 'running' | 'failed' | 'pending'>('all');

  const getStatusColor = (status: string) => {
    switch (status) {
      case 'completed': return 'bg-green-100 text-green-800';
      case 'running': return 'bg-blue-100 text-blue-800';
      case 'failed': return 'bg-red-100 text-red-800';
      case 'pending': return 'bg-yellow-100 text-yellow-800';
      default: return 'bg-gray-100 text-gray-800';
    }
  };

  const getStatusIcon = (status: string) => {
    switch (status) {
      case 'completed':
        return <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 12l2 2 4-4m6 2a9 9 0 11-18 0 9 9 0 0118 0z" />
        </svg>;
      case 'running':
        return <svg className="w-4 h-4 animate-spin" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15" />
        </svg>;
      case 'failed':
        return <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M10 14l2-2m0 0l2-2m-2 2l-2-2m2 2l2 2m7-2a9 9 0 11-18 0 9 9 0 0118 0z" />
        </svg>;
      case 'pending':
        return <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 8v4l3 3m6-3a9 9 0 11-18 0 9 9 0 0118 0z" />
        </svg>;
      default:
        return null;
    }
  };

  const filteredWorkflows = filter === 'all' 
    ? workflows 
    : workflows.filter(workflow => workflow.status === filter);

  const formatDate = (dateString: string) => {
    return new Date(dateString).toLocaleDateString('en-US', {
      year: 'numeric',
      month: 'short',
      day: 'numeric',
      hour: '2-digit',
      minute: '2-digit'
    });
  };

  return (
    <div className="min-h-screen bg-gray-50">
      <Navigation />
      
      <div className="container mx-auto px-4 py-8">
        <div className="mb-8">
          <h1 className="text-3xl font-bold text-gray-900 mb-4">Workflow Dashboard</h1>
          <p className="text-gray-600">Monitor and manage your bioinformatics workflows</p>
        </div>

        <div className="bg-white rounded-lg shadow-sm border border-gray-200 mb-6">
          <div className="p-6 border-b border-gray-200">
            <h2 className="text-lg font-semibold mb-4">Filter Workflows</h2>
            <div className="flex flex-wrap gap-2">
              {(['all', 'completed', 'running', 'failed', 'pending'] as const).map((status) => (
                <button
                  key={status}
                  onClick={() => setFilter(status)}
                  className={`px-4 py-2 rounded-lg text-sm font-medium transition-colors ${
                    filter === status
                      ? 'bg-blue-600 text-white'
                      : 'bg-gray-100 text-gray-700 hover:bg-gray-200'
                  }`}
                >
                  {status === 'all' ? 'All' : status.charAt(0).toUpperCase() + status.slice(1)}
                  {status !== 'all' && (
                    <span className="ml-2 text-xs">
                      ({workflows.filter(w => w.status === status).length})
                    </span>
                  )}
                </button>
              ))}
            </div>
          </div>

          <div className="p-6">
            {filteredWorkflows.length === 0 ? (
              <div className="text-center py-12">
                <svg className="mx-auto h-12 w-12 text-gray-400" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5H7a2 2 0 00-2 2v10a2 2 0 002 2h8a2 2 0 002-2V7a2 2 0 00-2-2h-2M9 5a2 2 0 002 2h2a2 2 0 002-2M9 5a2 2 0 012-2h2a2 2 0 012 2" />
                </svg>
                <h3 className="mt-2 text-sm font-medium text-gray-900">No workflows found</h3>
                <p className="mt-1 text-sm text-gray-500">
                  {filter === 'all' ? 'Create your first workflow to get started.' : `No workflows with ${filter} status.`}
                </p>
              </div>
            ) : (
              <div className="grid gap-4 md:grid-cols-2 lg:grid-cols-3">
                {filteredWorkflows.map((workflow) => (
                  <div
                    key={workflow.id}
                    className="border border-gray-200 rounded-lg p-4 hover:shadow-md transition-shadow cursor-pointer"
                  >
                    <div className="flex items-start justify-between mb-3">
                      <h3 className="text-lg font-semibold text-gray-900 truncate pr-2">
                        {workflow.name}
                      </h3>
                      <div className={`inline-flex items-center px-2.5 py-0.5 rounded-full text-xs font-medium ${getStatusColor(workflow.status)}`}>
                        {getStatusIcon(workflow.status)}
                        <span className="ml-1">{workflow.status}</span>
                      </div>
                    </div>
                    
                    <p className="text-sm text-gray-600 mb-3 line-clamp-2">
                      {workflow.description}
                    </p>
                    
                    <div className="flex items-center justify-between text-xs text-gray-500">
                      <span>{workflow.nodes_count} nodes</span>
                      <span>{formatDate(workflow.created_at)}</span>
                    </div>
                    
                    <div className="mt-4 flex space-x-2">
                      <button className="flex-1 bg-blue-600 text-white px-3 py-2 rounded text-sm hover:bg-blue-700 transition-colors">
                        View Details
                      </button>
                      {workflow.status === 'failed' && (
                        <button className="flex-1 bg-green-600 text-white px-3 py-2 rounded text-sm hover:bg-green-700 transition-colors">
                          Retry
                        </button>
                      )}
                      {workflow.status === 'running' && (
                        <button className="flex-1 bg-red-600 text-white px-3 py-2 rounded text-sm hover:bg-red-700 transition-colors">
                          Stop
                        </button>
                      )}
                    </div>
                  </div>
                ))}
              </div>
            )}
          </div>
        </div>

        <div className="bg-white rounded-lg shadow-sm border border-gray-200">
          <div className="p-6">
            <h2 className="text-lg font-semibold mb-4">Quick Stats</h2>
            <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
              <div className="text-center">
                <div className="text-2xl font-bold text-green-600">
                  {workflows.filter(w => w.status === 'completed').length}
                </div>
                <div className="text-sm text-gray-500">Completed</div>
              </div>
              <div className="text-center">
                <div className="text-2xl font-bold text-blue-600">
                  {workflows.filter(w => w.status === 'running').length}
                </div>
                <div className="text-sm text-gray-500">Running</div>
              </div>
              <div className="text-center">
                <div className="text-2xl font-bold text-red-600">
                  {workflows.filter(w => w.status === 'failed').length}
                </div>
                <div className="text-sm text-gray-500">Failed</div>
              </div>
              <div className="text-center">
                <div className="text-2xl font-bold text-yellow-600">
                  {workflows.filter(w => w.status === 'pending').length}
                </div>
                <div className="text-sm text-gray-500">Pending</div>
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}