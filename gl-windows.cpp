#include <assert.h>
#include <windows.h>
#include "glad.h"
#include "glad_wgl.h"

#include <iostream>
#include <string>
#include <thread>

static int width = 500;
static int height = 500;

static HWND windowHandle;
static HDC gdiDeviceContext;
static HGLRC glRenderContext;

static GLuint vao = 0;

static GLuint vertexShader = 0;
static GLuint fragmentShader = 0;
static GLuint shaderProgram = 0;

static std::string vertexShaderSource =
	"#version 430\n\
	layout(location = 0) in vec2 pos;\
	out vec2 p;\
	void main() {\
		gl_Position = vec4(pos, 0, 1);\
		p = pos * .5 + .5;\
	}";

static const std::string glslVersionString = "#version 430\n";
static std::string glslUniformStrings = "";
static std::string fragmentShaderSource =
"in vec2 p;\n\
out vec4 color;\n\
int i; float f;\n\
void main() {\n\
  color = vec4(p, 0, 1);\n";
static std::string fragmentShaderSourceTmp;

static bool running = true;
static bool inputThreadFlag = false;

/**
*
*/
static bool compileGLShader(std::string vertSrc, std::string fragSrc, GLuint* outVertexShader, GLuint* outFragmentShader) {
	const auto vertStr = vertSrc.c_str();
	const auto fragStr = fragSrc.c_str();
	const int vertStrLen = vertSrc.length();
	const int fragStrLen = fragSrc.length();
	GLuint v = glCreateShader(GL_VERTEX_SHADER);
	GLuint f = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(v, 1, &vertStr, &vertStrLen);
	glShaderSource(f, 1, &fragStr, &fragStrLen);
	glCompileShader(f);
	glCompileShader(v);
	GLint success;
	glGetShaderiv(f, GL_COMPILE_STATUS, &success);
	if (success != GL_TRUE) {
		GLchar errorString[1024];
		glGetShaderInfoLog(f, 1024, 0, errorString);
		std::cout << std::endl << errorString;
		return false;
	}
	*outVertexShader = v;
	*outFragmentShader = f;
	return true;
}

/**
*
*/
static void runGLShader() {
	glBindVertexArray(vao);
	glUseProgram(shaderProgram);
	glDisable(GL_DEPTH_TEST);
	glViewport(0, 0, width, height);
	glDrawArrays(GL_TRIANGLES, 0, 6);
}

/**
*
*/
static LRESULT CALLBACK onWindowMessage(
	HWND windowHandle,
	UINT uMsg,        // message identifier
	WPARAM wParam,    // first message parameter
	LPARAM lParam)    // second message parameter
{ 
	// https://docs.microsoft.com/en-us/windows/desktop/winmsg/about-messages-and-message-queues
	switch (uMsg) 
	{ 
		case WM_CREATE: 
			// Initialize the window. 
			return 0; 

		case WM_PAINT: 
			// Paint the window's client area. 
			runGLShader();
			SwapBuffers(GetDC(windowHandle));
			return 0; 

		case WM_SIZE: 
			// Set the size and position of the window. 
			return 0; 

		case WM_CHAR:
			if (wParam == VK_ESCAPE) {
				PostQuitMessage(0);
			}
			return 0;

		case WM_DESTROY: 
			// Clean up window-specific data objects. 
			PostQuitMessage(0); // puts a WM_QUIT message on the message queue, causing the message loop to end.
			return 0; 

		// 
		// Process other messages. 
		// 

		default: 
			return DefWindowProc(windowHandle, uMsg, wParam, lParam); 
	} 
	return 0; 
}

/**
*
*/
static int runWindowsMessageLoop() {
	// https://en.wikipedia.org/wiki/Message_loop_in_Microsoft_Windows
	MSG msg;
	BOOL bRet;
	while(running)
	{

		if (inputThreadFlag) {
			GLuint newVertexShader, newFragmentShader;
			if (compileGLShader(vertexShaderSource, glslVersionString + glslUniformStrings + fragmentShaderSourceTmp + "}", 
				&newVertexShader, &newFragmentShader)) 
			{
				fragmentShaderSource = fragmentShaderSourceTmp;
				glDeleteShader(fragmentShader);
				glDeleteShader(vertexShader);
				glDeleteProgram(shaderProgram);
				vertexShader = newVertexShader;
				fragmentShader = newFragmentShader;
				shaderProgram = glCreateProgram();
				glAttachShader(shaderProgram, fragmentShader);
				glAttachShader(shaderProgram, vertexShader);
				glLinkProgram(shaderProgram);
				std::cout << " [OK]" << std::endl;
			} else {
				std::cout << "  [CONTINUE AFTER ERROR]" << std::endl;
			}
			std::cout << "  " << std::flush; // indentation
			inputThreadFlag = false;
		}

		bRet = GetMessage(&msg, NULL, 0, 0);

		if (bRet > 0)  // (bRet > 0 indicates a message that must be processed.)
		{
			TranslateMessage(&msg); // Translates virtual-key messages into character messages. 
				// The character messages are posted to the calling thread's message queue, 
				// to be read the next time the thread calls the GetMessage or PeekMessage function.
			DispatchMessage(&msg); // Dispatches a message to a window procedure
		}
		else if (bRet < 0)  // (bRet == -1 indicates an error.)
		{
			// Handle or log the error; possibly exit.
			// ...
		}
		else  // (bRet == 0 indicates "exit program".)
		{
			break;
		}
	}
	return msg.wParam; // exit code from PostQuitMessage
}

/**
*
*/
static void createWindowsWindow(const char* title, int w, int h, HWND* outWindowHandle) {
	// Handle of exe or dll, depending on where call happens
	HINSTANCE moduleHandle = (HINSTANCE)GetModuleHandle(NULL);

	WNDCLASSEX wcx;
	// Fill in the window class structure with parameters 
	// that describe the main window. 
	wcx.cbSize = sizeof(wcx);          // size of structure 
	wcx.style = CS_OWNDC;              // Allocates a unique device context (DC) for each window in the class. 
		// https://docs.microsoft.com/en-us/windows/desktop/gdi/about-device-contexts
	wcx.lpfnWndProc = onWindowMessage; // points to window procedure 
	wcx.cbClsExtra = 0;                // no extra class memory 
	wcx.cbWndExtra = 0;                // no extra window memory 
	wcx.hInstance = moduleHandle;      // handle to instance 
	wcx.hIcon = LoadIcon(NULL, 
		IDI_APPLICATION);              // predefined app. icon 
	wcx.hCursor = LoadCursor(NULL, 
		IDC_ARROW);                    // predefined arrow 
	wcx.hbrBackground = static_cast<HBRUSH>(GetStockObject( 
		WHITE_BRUSH));                  // white background brush 
	wcx.lpszMenuName =  "MainMenu";    // name of menu resource 
	wcx.lpszClassName = "MainWClass";  // name of window class 
	wcx.hIconSm = static_cast<HICON>(LoadImage(moduleHandle, // small class icon 
		MAKEINTRESOURCE(5),
		IMAGE_ICON, 
		GetSystemMetrics(SM_CXSMICON), 
		GetSystemMetrics(SM_CYSMICON), 
		LR_DEFAULTCOLOR)); 

	// Register the window class. 
	RegisterClassEx(&wcx);

	// Now createWindow can be called
	HWND windowHandle = CreateWindow( 
		"MainWClass",        // name of window class 
		title,               // title-bar string 
		WS_OVERLAPPEDWINDOW, // top-level window 
		CW_USEDEFAULT,       // default horizontal position 
		CW_USEDEFAULT,       // default vertical position 
		w,                   // width 
		h,                   // height 
		(HWND) NULL,         // no owner window 
		(HMENU) NULL,        // use class menu 
		moduleHandle,        // handle to application instance 
		(LPVOID) NULL);      // no window-creation data
	*outWindowHandle = windowHandle;
}

/**
*
*/
typedef HGLRC (__stdcall * PFNWGLCREATECONTEXTATTRIBSARB) (HDC hDC, HGLRC hShareContext, const int *attribList);
typedef BOOL (__stdcall * PFNWGLCHOOSEPIXELFORMATARBPROC)(HDC hdc, const int * piAttribIList, const FLOAT * pfAttribFList, UINT nMaxFormats, int * piFormats, UINT * nNumFormats);

/**
*
*/
static HWND createDummyWindow() {
	DWORD windowExStyle, windowStyle;

	WNDCLASS	wc;						// Windows Class Structure
	RECT		WindowRect;				// Grabs Rectangle Upper Left / Lower Right Values
	WindowRect.left = 0L;
	WindowRect.right = 640L;
	WindowRect.top = 0L;
	WindowRect.bottom = 480L;

	HINSTANCE instance	= ::GetModuleHandle( NULL );				// Grab An Instance For Our Window
	wc.style			= CS_HREDRAW | CS_VREDRAW | CS_OWNDC;	// Redraw On Size, And Own DC For Window.
	wc.lpfnWndProc		= DefWindowProc;						// WndProc Handles Messages
	wc.cbClsExtra		= 0;									// No Extra Window Data
	wc.cbWndExtra		= 0;									// No Extra Window Data
	wc.hInstance		= instance;
	wc.hIcon			= ::LoadIcon( NULL, IDI_WINLOGO );		// Load The Default Icon
	wc.hCursor			= ::LoadCursor( NULL, IDC_ARROW );		// Load The Arrow Pointer
	wc.hbrBackground	= NULL;									// No Background Required For GL
	wc.lpszMenuName		= NULL;									// We Don't Want A Menu
	wc.lpszClassName	= TEXT("FLINTTEMP");

	if( ! ::RegisterClass( &wc ) ) {											// Attempt To Register The Window Class
		DWORD err = ::GetLastError();
		return 0;											
	}
	windowExStyle = WS_EX_APPWINDOW | WS_EX_WINDOWEDGE | WS_EX_ACCEPTFILES;		// Window Extended Style
	windowStyle = ( WS_OVERLAPPEDWINDOW & ~WS_THICKFRAME );

	::AdjustWindowRectEx( &WindowRect, windowStyle, FALSE, windowExStyle );

	return ::CreateWindowEx( windowExStyle, TEXT("FLINTTEMP"), TEXT("FLINT"), windowStyle, 0, 0, WindowRect.right-WindowRect.left, WindowRect.bottom-WindowRect.top, NULL, NULL, instance, 0 );
}

/**
*
*/
static bool getWglFunctionPointers( PFNWGLCREATECONTEXTATTRIBSARB *resultCreateContextAttribsFnPtr, PFNWGLCHOOSEPIXELFORMATARBPROC *resultChoosePixelFormatFnPtr ) {
	static PFNWGLCREATECONTEXTATTRIBSARB cachedCreateContextAttribsFnPtr = nullptr;
	static PFNWGLCHOOSEPIXELFORMATARBPROC cachedChoosePixelFormatFnPtr = nullptr;
	if( ! cachedCreateContextAttribsFnPtr || ! cachedChoosePixelFormatFnPtr ) {
		static PIXELFORMATDESCRIPTOR pfd = {
			sizeof(PIXELFORMATDESCRIPTOR),				// Size Of This Pixel Format Descriptor
			1,											// Version Number
			PFD_DRAW_TO_WINDOW |						// Format Must Support Window
			PFD_SUPPORT_OPENGL |						// Format Must Support OpenGL
			PFD_DOUBLEBUFFER,							// Must Support Double Buffering
			PFD_TYPE_RGBA,								// Request An RGBA Format
			32,											// Select Our Color Depth
			0, 0, 0, 0, 0, 0,							// Color Bits Ignored
			0,											// No Alpha Buffer
			0,											// Shift Bit Ignored
			0,											// No Accumulation Buffer
			0, 0, 0, 0,									// Accumulation Bits Ignored
			16,											// depth bits
			0,											// stencil bits
			0,											// No Auxiliary Buffer
			PFD_MAIN_PLANE,								// Main Drawing Layer
			0,											// Reserved
			0, 0, 0										// Layer Masks Ignored
		};

		HWND tempWindow = createDummyWindow();
		HDC tempDc = ::GetDC( tempWindow );
		auto pixelFormat = ::ChoosePixelFormat( tempDc, &pfd );
		if( pixelFormat == 0 ) {
			::ReleaseDC( tempWindow, tempDc );
			::DestroyWindow( tempWindow );
			::UnregisterClass( TEXT("FLINTTEMP"), ::GetModuleHandle( NULL ) );
			return false;
		}
		::SetPixelFormat( tempDc, pixelFormat, &pfd );
		auto tempCtx = ::wglCreateContext( tempDc ); 
		::wglMakeCurrent( tempDc, tempCtx );

		cachedCreateContextAttribsFnPtr = (PFNWGLCREATECONTEXTATTRIBSARB) ::wglGetProcAddress( "wglCreateContextAttribsARB" );
		cachedChoosePixelFormatFnPtr = (PFNWGLCHOOSEPIXELFORMATARBPROC) ::wglGetProcAddress( "wglChoosePixelFormatARB" );
		*resultCreateContextAttribsFnPtr = cachedCreateContextAttribsFnPtr;
		*resultChoosePixelFormatFnPtr = cachedChoosePixelFormatFnPtr;
		::wglMakeCurrent( NULL, NULL );
		::wglDeleteContext( tempCtx );

		::ReleaseDC( tempWindow, tempDc );
		::DestroyWindow( tempWindow );
		::UnregisterClass( TEXT("FLINTTEMP"), ::GetModuleHandle( NULL ) );

		if( ! cachedCreateContextAttribsFnPtr || ! cachedChoosePixelFormatFnPtr ) {
			return false;
		}
		else
			return true;
	}
	else {
		*resultCreateContextAttribsFnPtr = cachedCreateContextAttribsFnPtr;
		*resultChoosePixelFormatFnPtr = cachedChoosePixelFormatFnPtr;
		return cachedCreateContextAttribsFnPtr && cachedChoosePixelFormatFnPtr;
	}
}

/**
*
*/
static void createGLContext(HDC deviceContext, HGLRC* outGLContext) {
	// https://www.khronos.org/opengl/wiki/Creating_an_OpenGL_Context_(WGL)
	PIXELFORMATDESCRIPTOR pfd = {
		sizeof(PIXELFORMATDESCRIPTOR),
		1,
		PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER,    // Flags
		PFD_TYPE_RGBA,        // The kind of framebuffer. RGBA or palette.
		32,                   // Colordepth of the framebuffer.
		0, 0, 0, 0, 0, 0,
		0,
		0,
		0,
		0, 0, 0, 0,
		24,                   // Number of bits for the depthbuffer
		8,                    // Number of bits for the stencilbuffer
		0,                    // Number of Aux buffers in the framebuffer.
		PFD_MAIN_PLANE,
		0,
		0, 0, 0
	};
	int pixelFormat = ChoosePixelFormat(deviceContext, &pfd);
	SetPixelFormat(deviceContext, pixelFormat, &pfd);

	// https://www.khronos.org/registry/OpenGL/extensions/ARB/WGL_ARB_create_context.txt
	PFNWGLCREATECONTEXTATTRIBSARB wglCreateContextAttribsARBPtr = NULL;
	PFNWGLCHOOSEPIXELFORMATARBPROC wglChoosePixelFormatARBPtr = NULL;
	assert( getWglFunctionPointers( &wglCreateContextAttribsARBPtr, &wglChoosePixelFormatARBPtr ) );
	int attribList[] = {
		WGL_CONTEXT_MAJOR_VERSION_ARB, 4,
		WGL_CONTEXT_MINOR_VERSION_ARB, 0,
		WGL_CONTEXT_FLAGS_ARB, 0,
		WGL_CONTEXT_PROFILE_MASK_ARB, WGL_CONTEXT_CORE_PROFILE_BIT_ARB,
		0, 0
	};
	HGLRC glContext = wglCreateContextAttribsARBPtr(deviceContext, 0, attribList);

	wglMakeCurrent(deviceContext, glContext);
	*outGLContext = glContext;
}

/**
*
*/
static void loadGLFunctions() {
	gladLoadGL();
	// wglGetProcAddress();
}

/**
*
*/
static void createGLQuadVAO(GLuint* outVAO) {
	float quad[] = {
		// (x, y)
		-1.0f,  1.0f, // top left
		1.0f, -1.0f, // bottom right
		-1.0f, -1.0f, // bottom left
		-1.0f,  1.0f, // top left
		1.0f,  1.0f, // top right
		1.0f, -1.0f // bottom right
	};
	GLuint id;
	glGenVertexArrays(1, &id);
	glBindVertexArray(id);
	GLuint buffer;
	glGenBuffers(1, &buffer);
	glBindBuffer(GL_ARRAY_BUFFER, buffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(quad), quad, GL_STATIC_DRAW);
	for (int i = 0; i < 6; i++) {
		glVertexAttribPointer(i, 2, GL_FLOAT, GL_FALSE, 0, 0);
		glEnableVertexAttribArray(i);
	}
	*outVAO = id;
}

/**
*
*/
static void runReadPrintLoop() {
	std::cout << " ____________________________________________________________" << std::endl;
	std::cout << "|::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::|" << std::endl;
	std::cout << "|:::::::::::::::::: FRAGMENT SHADER EDITOR ::::::::::::::::::|" << std::endl;
	std::cout << "|::::::: works only in Microsoft shell, ESC to close ::::::::|" << std::endl;
	std::cout << std::endl;
	std::cout << glslVersionString + glslUniformStrings + fragmentShaderSource;
	const int pollInterval = 200; // milliseconds
	std::cout << "  " << std::flush; // indentation
	std::string in;
	while (running) {
		// poll standard input before reading so that it is non-blocking
		// this will not work for other shells than Microsoft
		if (WaitForSingleObject(GetStdHandle(STD_INPUT_HANDLE), pollInterval) == WAIT_OBJECT_0) {
			INPUT_RECORD buf[128]; DWORD readCount;
			if (ReadConsoleInput(GetStdHandle(STD_INPUT_HANDLE), buf, 128, &readCount)) {
				for (int i = 0; i < readCount; i++)  {
		            if (buf[i].EventType == KEY_EVENT) {
		            	if (buf[i].Event.KeyEvent.bKeyDown) {
		                    char c = buf[i].Event.KeyEvent.uChar.AsciiChar;
		                    if (c == '\r' || c == '\n') { // enter
		                    	if (in.size() > 0) {
									fragmentShaderSourceTmp = fragmentShaderSource + in;
									inputThreadFlag = true;
									in = "";
								}
		                    } else if (c == 8) { // backspace
		                    	if (in.size() > 0) {
		                    		std::string whitespace(in.size(), ' ');
		                    		std::cout << '\r' << "  " << whitespace << std::flush;
		                    		in.resize(in.size() - 1);
		                    		std::cout << '\r' << "  " << in << std::flush;
		                    	}
		                    } else if (c >= 32 && c <= 126) { // printable
		                    	in += c;
		                    	std::cout << c;
		                    }
		                    if (buf[i].Event.KeyEvent.wVirtualKeyCode == 27) { // escape
		                    	running = false;
		                    	break;
		                    }
		                }
		            }
		        }
			}
		}
	}
	std::cout << std::endl << "}";
}


// public:

/**
*
*/
void createGLContexts(void* outDeviceContext = 0, void* outRenderContext = 0) {
	createWindowsWindow("Shader Output", width, height, &windowHandle);
	gdiDeviceContext = GetDC(windowHandle);
	createGLContext(gdiDeviceContext, &glRenderContext);
	if (outDeviceContext && outRenderContext) {
		*((HDC*)outDeviceContext) = gdiDeviceContext;
		*((HGLRC*)outRenderContext) = glRenderContext;
	}

	loadGLFunctions();

	createGLQuadVAO(&vao);
}

/**
*
*/
void createGLBuffer(size_t bytes, void* outBuffer) {
	GLint maxSize; glGetIntegerv(GL_MAX_UNIFORM_BLOCK_SIZE, &maxSize);
	if (bytes > maxSize) {
		std::cout << bytes << " exceeds max size of " << maxSize << " for GL uniform buffer" << std::endl;
		return;
	}
	static GLuint bindingPoint = 0; // assignemt happens once at program start, then value persists

	std::string bindingPointString = std::to_string(bindingPoint);
	std::string uniformName = "Buf" + bindingPointString;
	size_t words = bytes / 4;

	// to access each scalar in array, use std140 memory layout rules and access through vec4
	// https://www.khronos.org/opengl/wiki/Interface_Block_(GLSL)#Memory_layout
	glslUniformStrings += "layout(std140, binding=" + bindingPointString + ") uniform " + uniformName
		+ " {\n  uvec4 buf" + bindingPointString + "[" + std::to_string(words/4) + "];"
		+ " // " + std::to_string(words) + " uints, " + std::to_string(bytes) + " bytes"
		+ "\n};\n";
	fragmentShaderSource = fragmentShaderSource
		+ "  i = int(floor(p.x * " + std::to_string(words) + "));\n"
		+ "  f = float(buf" + bindingPointString + "[i/4][i%4]) / 4294967296.;\n"
		+ "  if(p.y < f) color=vec4(1); else color=vec4(0);\n";
	std::string shaderSource = glslVersionString + glslUniformStrings + fragmentShaderSource;
	compileGLShader(vertexShaderSource, shaderSource + "}", &vertexShader, &fragmentShader);
	shaderProgram = glCreateProgram();
	glAttachShader(shaderProgram, fragmentShader);
	glAttachShader(shaderProgram, vertexShader);
	glLinkProgram(shaderProgram);

	GLuint ubo = 0;
	glGenBuffers(1, &ubo);
	glBindBuffer(GL_UNIFORM_BUFFER, ubo);
	glBufferData(GL_UNIFORM_BUFFER, bytes, NULL, GL_DYNAMIC_DRAW);
	glBindBuffer(GL_UNIFORM_BUFFER, 0);
	GLuint index = glGetUniformBlockIndex(shaderProgram, uniformName.c_str());
	if (index != GL_INVALID_INDEX) {
		glUniformBlockBinding(shaderProgram, index, bindingPoint);
		glBindBufferBase(GL_UNIFORM_BUFFER, bindingPoint, ubo);
	}
	bindingPoint++;
	*((GLuint*)outBuffer) = ubo;
	//TODO buffer cleanup at end of render loop
}

/**
*
*/
void createGLImage(size_t w, size_t h, void* outImageHandle, void* data = 0) {
	static int i = 0;
	glslUniformStrings += "layout(location = 0) uniform sampler2D img" + std::to_string(i) + ";\n";
	fragmentShaderSource += "  color = texture(img" + std::to_string(i) + ", p);\n";
	std::string shaderSource = glslVersionString + glslUniformStrings + fragmentShaderSource;
	compileGLShader(vertexShaderSource, shaderSource + "}", &vertexShader, &fragmentShader);
	shaderProgram = glCreateProgram();
	glAttachShader(shaderProgram, fragmentShader);
	glAttachShader(shaderProgram, vertexShader);
	glLinkProgram(shaderProgram);

	i++;

	GLuint tex;
	glGenTextures(1, &tex),
	glBindTexture(GL_TEXTURE_2D, tex);

	// Use GL_R32_F as internal format (i.e. GL_RED as general format) because only greyscale needed
	glTexImage2D(GL_TEXTURE_2D, 0/*mip*/, GL_R32F, w, h, 0, GL_RED, GL_FLOAT, (GLvoid*)data);
	*((GLuint*)outImageHandle) = tex;
}

/**
*
*/
void runGLRenderLoop() {
	ShowWindow(windowHandle, SW_SHOWNORMAL);
	UpdateWindow(windowHandle);

	compileGLShader(vertexShaderSource, glslVersionString + glslUniformStrings + fragmentShaderSource + "}",
		&vertexShader, &fragmentShader);
	shaderProgram = glCreateProgram();
	glAttachShader(shaderProgram, fragmentShader);
	glAttachShader(shaderProgram, vertexShader);
	glLinkProgram(shaderProgram);

	std::thread inputThread(runReadPrintLoop);

	runWindowsMessageLoop();

	// window closed now, wait for input thread before cleanup
	running = false; inputThread.join();

	glDeleteShader(fragmentShader);
	glDeleteShader(vertexShader);
	glDeleteProgram(shaderProgram);

	glDeleteVertexArrays(1, &vao);

	// make the rendering context not current before deleting it
	wglMakeCurrent(NULL, NULL);
	wglDeleteContext(glRenderContext);
}



// https://msdn.microsoft.com/de-de/library/windows/desktop/ms633504(v=vs.85).aspx
// GetDesktopWindow 
// Retrieves a handle to the desktop window. The desktop window covers the entire screen. 
// The desktop window is the area on top of which other windows are painted.