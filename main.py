# FastAPI main application - Pure API Backend
import os
from fastapi import FastAPI, Request, UploadFile, File, Form, HTTPException, Depends
from fastapi.responses import JSONResponse, FileResponse
from fastapi.middleware.cors import CORSMiddleware
from dotenv import load_dotenv
from werkzeug.utils import secure_filename

# Import helper and routers
import helper

# Load environment variables
load_dotenv()

# Create FastAPI app
app = FastAPI(
    title="CCST - Cloud Control and Simulation Toolbox",
    description="Web-based control systems simulation and aerospace engineering toolkit",
    version="2.0.0"
)

# CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost:3000",  # Next.js dev server
        "http://127.0.0.1:3000",  # Alternative localhost
        "*"  # Allow all for development (configure for production)
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Configuration
APP_ROOT = os.path.dirname(os.path.abspath(__file__))
UPLOAD_FOLDER = os.path.join(APP_ROOT, 'uploads')
ALLOWED_EXTENSIONS = {'png', 'jpg', 'jpeg', 'gif', 'bmp'}
MAX_FILE_SIZE = 16 * 1024 * 1024  # 16MB
API_KEY = os.getenv('API_KEY', None)

# Ensure upload folder exists
os.makedirs(UPLOAD_FOLDER, exist_ok=True)

# Initialize states
states = {}
helper.save_states(states)


# Security: API Key dependency
async def verify_api_key(request: Request):
    """Verify API key if one is configured"""
    if API_KEY:
        provided_key = request.headers.get('X-API-Key') or request.query_params.get('api_key')
        if not provided_key or provided_key != API_KEY:
            raise HTTPException(status_code=401, detail="Invalid or missing API key")
    return True


def allowed_file(filename: str) -> bool:
    """Check if file extension is allowed"""
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


# API Routes
@app.get("/")
async def root():
    """API root endpoint"""
    return {
        "message": "CCST API Backend",
        "version": "2.0.0",
        "docs": "/docs",
        "frontend": "http://localhost:3000 (Next.js)"
    }


@app.post("/upload")
async def upload_file(
    request: Request,
    file: UploadFile = File(...),
    inputTitle: str = Form(...),
    _: bool = Depends(verify_api_key)
):
    """Upload image file"""
    if not file.filename:
        raise HTTPException(status_code=400, detail="No file selected")

    if not allowed_file(file.filename):
        raise HTTPException(
            status_code=400,
            detail=f"File type not allowed. Allowed types: {', '.join(ALLOWED_EXTENSIONS)}"
        )

    # Secure the filename
    filename = secure_filename(file.filename)
    file_path = os.path.join(UPLOAD_FOLDER, filename)

    # Save file
    try:
        with open(file_path, "wb") as f:
            content = await file.read()
            f.write(content)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error saving file: {str(e)}")

    # Save state
    states = helper.load_states()
    states[inputTitle] = {
        'value': filename,
        'meta': {'what': 'image', 'value': filename}
    }
    helper.save_states(states)

    return {
        "message": "File uploaded successfully",
        "filename": filename,
        "title": inputTitle
    }


@app.get("/uploads/{filename}")
async def send_file(filename: str):
    """Serve uploaded files"""
    file_path = os.path.join(UPLOAD_FOLDER, filename)
    if not os.path.exists(file_path):
        raise HTTPException(status_code=404, detail="File not found")
    return FileResponse(file_path)


@app.get("/gdisplay")
async def gdisplay(name: str):
    """Retrieve variable metadata for graphical display"""
    states = helper.load_states()
    if name not in states:
        raise HTTPException(status_code=404, detail=f"Variable '{name}' not found")
    return states[name]['meta']


@app.get("/api/ping")
async def ping():
    """Simple ping endpoint for testing"""
    return {"message": "pong"}


@app.get("/api/display")
async def display(name: str):
    """Display variable value"""
    states = helper.load_states()
    if name not in states:
        raise HTTPException(status_code=404, detail=f"Variable '{name}' not found")

    meta = states[name]['meta']
    if meta['what'] == 'scalar':
        return {"value": states[name]['value']}
    elif meta['what'] == 'image':
        return {"message": f"To show the image type: gdisplay {name}", "value": meta['value']}
    elif meta['what'] == 'matrix':
        return {"value": states[name]['value'].tolist()}
    else:
        return {"value": states[name]['value'].tolist()}


# Error handlers
@app.exception_handler(404)
async def not_found(request: Request, exc: HTTPException):
    """Custom 404 handler"""
    return JSONResponse(
        status_code=404,
        content={"detail": "Not found", "path": request.url.path}
    )


# Import and include routers (we'll create these next)
# Note: We'll import these after creating the router files
try:
    from control_router import router as control_router
    app.include_router(control_router, prefix="/api/control", tags=["control"])
except ImportError:
    print("Warning: control_router not found, skipping...")

try:
    from data_router import router as data_router
    app.include_router(data_router, prefix="/api/data", tags=["data"])
except ImportError:
    print("Warning: data_router not found, skipping...")

try:
    from gnc_router import router as gnc_router
    app.include_router(gnc_router, prefix="/api/gnc", tags=["gnc"])
except ImportError:
    print("Warning: gnc_router not found, skipping...")


# Note: Use FastAPI CLI to run the application:
#   Development: uv run fastapi dev main.py
#   Production:  uv run fastapi run main.py --host 0.0.0.0 --port $PORT
