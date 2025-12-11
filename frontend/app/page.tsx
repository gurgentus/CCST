'use client';

import { useState } from 'react';
import Terminal from '@/components/Terminal';
import MatrixOutput from '@/components/MatrixOutput';
import { DisplayData } from '@/types/terminal';

export default function Home() {
  const [displayData, setDisplayData] = useState<DisplayData | null>(null);

  const handleTerminalOutput = (output: string) => {
    try {
      const data = JSON.parse(output);
      if (data.what) {
        setDisplayData(data);
      }
    } catch (e) {
      // Not JSON, ignore
    }
  };

  return (
    <main className="min-h-screen bg-gradient-to-br from-gray-900 via-gray-800 to-gray-900">
      {/* Header */}
      <header className="bg-gray-800 border-b border-gray-700 shadow-lg">
        <div className="container mx-auto px-4 py-6">
          <h1 className="text-3xl font-bold bg-gradient-to-r from-cyan-400 to-blue-500 bg-clip-text text-transparent">
            Cloud Controls and Simulation Toolbox
          </h1>
          <p className="text-gray-400 mt-2">
            FastAPI Backend + Next.js Frontend
          </p>
        </div>
      </header>

      {/* Main Content */}
      <div className="container mx-auto px-4 py-8">
        {/* Terminal Section */}
        <section className="mb-8">
          <h2 className="text-2xl font-bold mb-4 text-gray-200">Terminal</h2>
          <Terminal onOutput={handleTerminalOutput} />
        </section>

        {/* Output Section */}
        <section>
          <MatrixOutput data={displayData} />
        </section>

        {/* Quick Start Guide */}
        <section className="mt-12 bg-gray-800 border border-gray-700 rounded-lg p-6">
          <h3 className="text-xl font-bold mb-4 text-cyan-400">Quick Start</h3>
          <div className="grid md:grid-cols-2 gap-6 text-sm">
            <div>
              <h4 className="font-semibold text-gray-300 mb-2">Basic Commands</h4>
              <ul className="space-y-1 text-gray-400 font-mono">
                <li>{'>'} ping</li>
                <li>{'>'} help</li>
                <li>{'>'} matrix A [[1,2],[3,4]]</li>
                <li>{'>'} gdisplay A</li>
              </ul>
            </div>
            <div>
              <h4 className="font-semibold text-gray-300 mb-2">Control Systems</h4>
              <ul className="space-y-1 text-gray-400 font-mono">
                <li>{'>'} modes A V E</li>
                <li>{'>'} ss A B C D sys1</li>
                <li>{'>'} feedback G K sys_fb</li>
              </ul>
            </div>
          </div>
        </section>
      </div>

      {/* Footer */}
      <footer className="bg-gray-800 border-t border-gray-700 mt-12">
        <div className="container mx-auto px-4 py-6 text-center text-gray-500 text-sm">
          <p>© 2025 Cloud Controls and Simulation Toolbox</p>
          <p className="mt-1">
            Built with Next.js 15, React 19, TypeScript & Tailwind CSS
          </p>
        </div>
      </footer>
    </main>
  );
}
